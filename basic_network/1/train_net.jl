using Main.MOD_aeif
using Statistics
using LinearAlgebra
using MAT
using Random
using SparseArrays
using HDF5
include("parameters.jl")

function train_net(rcd_name)

println("Setting up neurons and networks")

T = 8000
rcd_buffer_x = 1 #the length of every rcd for spikes
rcd_buffer_w = 100 #the length of every rcd for weights
dt = 1e-4
frames = Int(round(T/dt))
buffer_frames_x = Int(round(rcd_buffer_x/dt))
buffer_frames_w = Int(round(rcd_buffer_w/dt))
interval_x = 1#Int(round(1/dt))

E = MOD_aeif.NeuronArray(N = Ne)
E.param.μθ = μθe
E.param.σθ = σθe
E.param.V_reset = V_reset_e
E.pastfires = zeros(Ne, Int(round(1.5e-3/dt)))
E.param.rt = rt_e
E.param.ηIP = ηIP_e

I = MOD_aeif.NeuronArray(N = Ni)
I.param.μθ = μθi
I.param.σθ = σθi
I.param.V_reset = V_reset_i
I.pastfires = zeros(Ni, Int(round(1.5e-3/dt)))
I.param.rt = rt_i
I.param.ηIP = ηIP_i

EE = MOD_aeif.Network(pre = E, post = E)
EE.w = sprand(Bool, Ne, Ne, 1.0)
EE.sparse = sprand(Bool, Ne, Ne, p_ee)
EE.sparse = EE.sparse - Diagonal(EE.sparse)
while minimum(EE.sparse*ones(Ne, 1)) == 0
    EE.sparse = sprand(Bool, Ne, Ne, p_ee)
    EE.sparse = EE.sparse - Diagonal(EE.sparse)
end
EE.w = EE.w .* EE.sparse
EE.sum_post = max.(0, broadcast(+, μTee, σTee*randn(EE.post.N)))
EE.w = EE.w .* (EE.sum_post * ones(1, Ne)) ./ (EE.w * ones(Ne, Ne))
EE.w = EE.w .* EE.sparse
EE.wt = SparseMatrixCSC(EE.w')
EE.sum_pre = EE.w' * ones(EE.post.N)
EE.sum_post_fixed = copy(EE.sum_post)
EE.sum_pre_fixed = copy(EE.sum_pre)
EE.param.exci_inh = 1
EE.param.plastic = true
EE.param.plasticity_rule = "eSTDP"
EE.param.delay = Int(round(d_ee/dt))
EE.param.A1 = AeLTP
EE.param.A2 = AeLTD
EE.param.τ1 = τLTP
EE.param.τ2 = τLTD

EI = MOD_aeif.Network(pre = E, post = I)
EI.w = sprand(Bool, Ni, Ne, 1.0)
EI.sparse = sprand(Bool, Ni, Ne, p_ei)
while minimum(EI.sparse*ones(Ne, 1)) == 0
    EI.sparse = sprand(Bool, Ni, Ne, p_ei)
end
EI.w = EI.w .* EI.sparse
EI.sum_post = max.(0, broadcast(+, μTei, σTei*randn(EI.post.N)))
EI.w = EI.w .* (EI.sum_post * ones(1, Ne)) ./ (EI.w * ones(Ne, Ne))
EI.w = EI.w .* EI.sparse
EI.wt = SparseMatrixCSC(EI.w')
EI.sum_pre = EI.w' * ones(EI.post.N)
EI.sum_post_fixed = copy(EI.sum_post)
EI.sum_pre_fixed = copy(EI.sum_pre)
EI.param.exci_inh = 1
EI.param.plastic = false
EI.param.delay = Int(round(d_ei/dt))

IE = MOD_aeif.Network(pre = I, post = E)
IE.w = sprand(Bool, Ne, Ni, 1.0)
IE.sparse = sprand(Bool, Ne, Ni, p_ie)
while minimum(IE.sparse*ones(Ni, 1)) == 0
    IE.sparse = sprand(Bool, Ne, Ni, p_ie)
end
IE.w = IE.w .* IE.sparse
IE.sum_post = max.(0, broadcast(+, μTie, σTie*randn(IE.post.N)))
IE.w = IE.w .* (IE.sum_post * ones(1, Ni)) ./ (IE.w * ones(Ni, Ni))
IE.w = IE.w .* IE.sparse
IE.wt = SparseMatrixCSC(IE.w')
IE.sum_pre = IE.w' * ones(IE.post.N)
IE.sum_post_fixed = copy(IE.sum_post)
IE.sum_pre_fixed = copy(IE.sum_pre)
IE.param.exci_inh = -1
IE.param.plastic = true
IE.param.plasticity_rule = "iSTDP"
IE.param.delay = Int(round(d_ie/dt))
IE.param.A1 = Aipre
IE.param.A2 = Aipost
IE.param.τ1 = τpre
IE.param.τ2 = τpost
IE.param.LTDα = LTDalpha

II = MOD_aeif.Network(pre = I, post = I)
II.w = sprand(Bool, Ni, Ni, 1.0)
II.sparse = sprand(Bool, Ni, Ni, p_ii)
II.sparse = II.sparse - Diagonal(II.sparse)
while minimum(II.sparse*ones(Ni, 1)) == 0
    II.sparse = sprand(Bool, Ni, Ni, p_ii)
    II.sparse = II.sparse - Diagonal(II.sparse)
end
II.w = II.w .* II.sparse
II.sum_post = max.(0, broadcast(+, μTii, σTii*randn(II.post.N)))
II.w = II.w .* (II.sum_post * ones(1, Ni)) ./ (II.w * ones(Ni, Ni))
II.w = II.w .* II.sparse
II.wt = SparseMatrixCSC(II.w')
II.sum_pre = II.w' * ones(II.post.N)
II.sum_post_fixed = copy(II.sum_post)
II.sum_pre_fixed = copy(II.sum_pre)
II.param.exci_inh = -1
II.param.plastic = false
II.param.delay = Int(round(d_ii/dt))

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
