using Main.MOD_aeif
using Statistics
using LinearAlgebra
using MAT
using Random
using SparseArrays
using HDF5

include("parameters.jl")
include("retrieve_sp.jl")

function train_net(rcd_name)

println("Setting up neurons and networks")

T = 17100
rcd_buffer_x = 1 #the length of every rcd for spikes
rcd_buffer_w = 100 #the length of every rcd for weights
dt = 1e-4
frames = Int(round(T/dt))
buffer_frames_x = Int(round(rcd_buffer_x/dt))
buffer_frames_w = Int(round(rcd_buffer_w/dt))
interval_x = 1#Int(round(1/dt))

ckpt = 80
network_path = string("../../basic_network/1/", rcd_name, "/")

ET = h5read(string(network_path, "all_w.h5"), string(ckpt, "/ET"))
IT = h5read(string(network_path, "all_w.h5"), string(ckpt, "/IT"))
EI_w = retrieve_sp(h5read(string(network_path, "init_condition_w.h5"), "EI"))
II_w = retrieve_sp(h5read(string(network_path, "init_condition_w.h5"), "II"))

E = MOD_aeif.NeuronArray(N = Ne)
E.T = ET
E.param.V_reset = V_reset_e
E.pastfires = zeros(Ne, Int(round(1.5e-3/dt)))
E.param.rt = rt_e
E.param.ηIP = 0
E.param.corr_t = corr_t_e
E.param.gap_t = gap_t_e
E.param.ξext = ξext_e
E.param.τext = τext_e
E.param.time_offset = 0

I = MOD_aeif.NeuronArray(N = Ni)
I.T = IT
I.param.V_reset = V_reset_i
I.pastfires = zeros(Ni, Int(round(1.5e-3/dt)))
I.param.rt = rt_i
I.param.ηIP = 0
I.M = 0
I.param.ξext = ξext_i
I.param.τext = τext_i
I.param.time_offset = 0

EI = MOD_aeif.Network(pre = E, post = I)
EI.w = EI_w
EI.sparse = sign.(EI.w)
EI.wt = SparseMatrixCSC(EI.w')
EI.sum_post = EI.w * ones(EI.pre.N)
EI.sum_pre = EI.w' * ones(EI.post.N)
EI.sum_post_fixed = copy(EI.sum_post)
EI.sum_pre_fixed = copy(EI.sum_pre)
EI.param.exci_inh = 1
EI.param.plastic = false
EI.param.delay = Int(round(d_ei/dt))

II = MOD_aeif.Network(pre = I, post = I)
II.w = II_w
II.sparse = sign.(II.w)
II.wt = SparseMatrixCSC(II.w')
II.sum_post = II.w * ones(II.pre.N)
II.sum_pre = II.w' * ones(II.post.N)
II.sum_post_fixed = copy(II.sum_post)
II.sum_pre_fixed = copy(II.sum_pre)
II.param.exci_inh = -1
II.param.plastic = false
II.param.delay = Int(round(d_ii/dt))

EE = MOD_aeif.Network(pre = E, post = E)
if ckpt == 0
    EE.w = retrieve_sp(h5read(string(network_path, "init_condition_w.h5"), "EE"))
else
    EE.w = retrieve_sp(h5read(string(network_path, "all_w.h5"), string(ckpt, "/EE")))
end
EE.sparse = sign.(EE.w)
EE.wt = SparseMatrixCSC(EE.w')
EE.sum_post = EE.w * ones(EE.pre.N)
EE.sum_pre = EE.w' * ones(EE.post.N)
EE.sum_post_fixed = copy(EE.sum_post)
EE.sum_pre_fixed = copy(EE.sum_pre)
EE.param.exci_inh = 1
EE.param.plastic = false
EE.param.plasticity_rule = "eSTDP"
EE.param.delay = Int(round(d_ee/dt))
EE.param.A1 = AeLTP
EE.param.A2 = AeLTD
EE.param.τ1 = τLTP
EE.param.τ2 = τLTD

IE = MOD_aeif.Network(pre = I, post = E)
if ckpt == 0
    IE.w = retrieve_sp(h5read(string(network_path, "init_condition_w.h5"), "IE"))
else
    IE.w = retrieve_sp(h5read(string(network_path, "all_w.h5"), string(ckpt, "/IE")))
end
IE.sparse = sign.(IE.w)
IE.wt = SparseMatrixCSC(IE.w')
IE.sum_post = IE.w * ones(IE.pre.N)
IE.sum_pre = IE.w' * ones(IE.post.N)
IE.sum_post_fixed = copy(IE.sum_post)
IE.sum_pre_fixed = copy(IE.sum_pre)
IE.param.exci_inh = -1
IE.param.plastic = false
IE.param.plasticity_rule = "iSTDP"
IE.param.delay = Int(round(d_ie/dt))
IE.param.A1 = Aipre
IE.param.A2 = Aipost
IE.param.τ1 = τpre
IE.param.τ2 = τpost
IE.param.LTDα = LTDalpha

system = MOD_aeif.the_whole_system(dt = dt, E = E, I = I, EE = EE, EI = EI,
    IE = IE, II = II, buffer_frames_x = buffer_frames_x, buffer_frames_w = buffer_frames_w)

rcdx = MOD_aeif.Record_x(frames = buffer_frames_x, interval_x = interval_x, Ne = Ne, Ni = Ni)

rcdw = MOD_aeif.Record_w(frames = buffer_frames_w, Ne = Ne, Ni = Ni)

MOD_aeif.initialize_system!(system)

println("Building output system")

mkdir(string(rcd_name))

MOD_aeif.record_w_init_cond(string(rcd_name, "/init_condition_w.h5"), system)
h5open(string(rcd_name, "/all_x.h5"), "w")
h5open(string(rcd_name, "/all_w.h5"), "w")
h5open(string(rcd_name, "/rcd_net_stat.h5"), "w")

println("Start running simulation")

for frame = 1:frames
    if ceil(frame*dt) == floor(frame*dt)
        println(frame*dt, " seconds")
        if ceil(frame*dt) == 100
            println("Plasticity turned on")
            EE.param.plastic = true
            IE.param.plastic = true
            E.T = ET
            I.T = IT
            E.param.ηIP = ηIP_e
            I.param.ηIP = ηIP_i
            system.EE.param.plastic = true
            system.IE.param.plastic = true
            system.E.T = ET
            system.I.T = IT
            system.E.param.ηIP = ηIP_e
            system.I.param.ηIP = ηIP_i
        end
    end
    MOD_aeif.run_process!(system, rcdx, rcdw, frame)

    if ceil(frame/buffer_frames_x) == floor(frame/buffer_frames_x)
        segment_no_x = Int(round(frame/buffer_frames_x))

        rcd_std_EE = std(system.EE.w[system.EE.sparse .== 1])
        rcd_std_IE = std(system.IE.w[system.IE.sparse .== 1])
        rcd_max_EE = maximum(system.EE.w[system.EE.sparse .== 1])
        rcd_max_IE = maximum(system.IE.w[system.IE.sparse .== 1])
        rcd_rate_E = sum(rcdx.xe)/Ne/rcd_buffer_x
        rcd_rate_I = sum(rcdx.xi)/Ni/rcd_buffer_x
        rcd_max_ET = maximum(system.E.T)
        rcd_min_ET = minimum(system.E.T)
        rcd_max_IT = maximum(system.I.T)
        rcd_min_IT = minimum(system.I.T)

        h5write(string(rcd_name, "/rcd_net_stat.h5"), string(segment_no_x, "/std_EE"), rcd_std_EE)
        h5write(string(rcd_name, "/rcd_net_stat.h5"), string(segment_no_x, "/max_EE"), rcd_max_EE)
        h5write(string(rcd_name, "/rcd_net_stat.h5"), string(segment_no_x, "/rate_E"), rcd_rate_E)
        h5write(string(rcd_name, "/rcd_net_stat.h5"), string(segment_no_x, "/rate_I"), rcd_rate_I)
        h5write(string(rcd_name, "/rcd_net_stat.h5"), string(segment_no_x, "/max_ET"), rcd_max_ET)
        h5write(string(rcd_name, "/rcd_net_stat.h5"), string(segment_no_x, "/min_ET"), rcd_min_ET)
        h5write(string(rcd_name, "/rcd_net_stat.h5"), string(segment_no_x, "/max_IT"), rcd_max_IT)
        h5write(string(rcd_name, "/rcd_net_stat.h5"), string(segment_no_x, "/min_IT"), rcd_min_IT)
    end
    if ceil(frame/buffer_frames_w) == floor(frame/buffer_frames_w)
        segment_no_w = Int(round(frame/buffer_frames_w))
        MOD_aeif.save_rcdw_hdf5(string(rcd_name, "/all_w.h5"), segment_no_w, rcdw)
        MOD_aeif.initialize_rcdw!(rcdw)
    end
    if ceil(frame/buffer_frames_x) == floor(frame/buffer_frames_x)
        segment_no_x = Int(round(frame/buffer_frames_x))
        MOD_aeif.save_rcdx_hdf5(string(rcd_name, "/all_x.h5"), segment_no_x, rcdx)
        MOD_aeif.initialize_rcdx!(rcdx)
        MOD_aeif.take_record_x!(0, system, rcdx, buffer_frames_x)

        if rcd_std_EE > 20 || rcd_max_EE > 500 || rcd_rate_E > 100
            file_break = matopen(string(rcd_name, "/blow_up.mat"), "w")
            write(file_break, "std_EE", rcd_std_EE)
            write(file_break, "max_EE", rcd_max_EE)
            write(file_break, "rate_E", rcd_rate_E)
            write(file_break, "time", frame*dt)
            close(file_break)
            break
        end

    end
end

end
