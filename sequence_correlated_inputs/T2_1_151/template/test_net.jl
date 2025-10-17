using Main.MODT_aeif
using Statistics
using LinearAlgebra
using MAT
using Random
using SparseArrays
using HDF5

include("parameters.jl")
include("retrieve_sp.jl")

function test_net(rcd_name)

println("Setting up neurons and networks.")

T = 43
rcd_buffer_x = 0.1 #the length of every rcd for spikes
rcd_buffer_w = T #the length of every rcd for weights
dt = 1e-4
frames = Int(round(T/dt))
buffer_frames_x = Int(round(rcd_buffer_x/dt))
buffer_frames_w = Int(round(rcd_buffer_w/dt))
interval_x = 1#Int(round(1/dt))

network_path = "../../../correlated_inputs/2/1/"
ckpt = 151

ET = h5read(string(network_path, "all_w.h5"), string(ckpt, "/ET"))
IT = h5read(string(network_path, "all_w.h5"), string(ckpt, "/IT"))
EI_w = retrieve_sp(h5read(string(network_path, "init_condition_w.h5"), "EI"))
II_w = retrieve_sp(h5read(string(network_path, "init_condition_w.h5"), "II"))

E = MODT_aeif.NeuronArray(N = Ne)
E.T = ET
E.param.V_reset = V_reset_e

I = MODT_aeif.NeuronArray(N = Ni)
I.T = IT
I.param.V_reset = V_reset_i

EI = MODT_aeif.Network(pre = E, post = I)
EI.w = EI_w
EI.sparse = sign.(EI.w)
EI.wt = SparseMatrixCSC(EI.w')
EI.param.exci_inh = 1
EI.param.plastic = false
EI.param.delay = Int(round(d_ei/dt))

II = MODT_aeif.Network(pre = I, post = I)
II.w = II_w
II.sparse = sign.(II.w)
II.wt = SparseMatrixCSC(II.w')
II.param.exci_inh = -1
II.param.plastic = false
II.param.delay = Int(round(d_ii/dt))

EE = MODT_aeif.Network(pre = E, post = E)
if ckpt == 0
    EE.w = retrieve_sp(h5read(string(network_path, "init_condition_w.h5"), "EE"))
else
    EE.w = retrieve_sp(h5read(string(network_path, "all_w.h5"), string(ckpt, "/EE")))
end
EE.sparse = sign.(EE.w)
EE.wt = SparseMatrixCSC(EE.w')
EE.param.exci_inh = 1
EE.param.plastic = false
EE.param.plasticity_rule = "eSTDP"
EE.param.delay = Int(round(d_ee/dt))
try
    E.inject_point = ((1:Ne) .== parse(Int, pwd()[end-3:end]))
catch
    try
        E.inject_point = ((1:Ne) .== parse(Int, pwd()[end-2:end]))
    catch
        try
            E.inject_point = ((1:Ne) .== parse(Int, pwd()[end-1:end]))
        catch
            E.inject_point = ((1:Ne) .== parse(Int, pwd()[end]))
        end
    end
end

IE = MODT_aeif.Network(pre = I, post = E)
if ckpt == 0
    IE.w = retrieve_sp(h5read(string(network_path, "init_condition_w.h5"), "IE"))
else
    IE.w = retrieve_sp(h5read(string(network_path, "all_w.h5"), string(ckpt, "/IE")))
end
IE.sparse = sign.(IE.w)
IE.wt = SparseMatrixCSC(IE.w')
IE.param.exci_inh = -1
IE.param.plastic = false
IE.param.plasticity_rule = "iSTDP"
IE.param.delay = Int(round(d_ie/dt))

println("Checkpoint ", ckpt, ", allocating memory for records")

system = MODT_aeif.the_whole_system(dt = dt, E = E, I = I, EE = EE, EI = EI,
    IE = IE, II = II, buffer_frames_x = buffer_frames_x, buffer_frames_w = buffer_frames_w)

rcdx = MODT_aeif.Record_x(frames = buffer_frames_x, interval_x = interval_x, Ne = Ne, Ni = Ni)

rcdw = MODT_aeif.Record_w(frames = buffer_frames_w, Ne = Ne, Ni = Ni)

MODT_aeif.initialize_test_system!(system)

println("Building output system")

mkdir(string(rcd_name))

MODT_aeif.record_w_init_cond(string(rcd_name, "/init_condition_w.h5"), system)
h5open(string(rcd_name, "/all_x.h5"), "w")
h5open(string(rcd_name, "/all_w.h5"), "w")

println("Checkpoint ", ckpt, ", start running simulation")

for frame = 1:frames
    if ceil(frame*dt) == floor(frame*dt)
        println("Checkpoint ", ckpt, ", ", frame*dt, " seconds")
    end
    MODT_aeif.run_process!(system, rcdx, rcdw, frame)

    if ceil(frame/buffer_frames_w) == floor(frame/buffer_frames_w)
        segment_no_w = Int(round(frame/buffer_frames_w))
        MODT_aeif.save_rcdw_hdf5(string(rcd_name, "/all_w.h5"), segment_no_w, rcdw)
        MODT_aeif.initialize_rcdw!(rcdw)
    end
    if ceil(frame/buffer_frames_x) == floor(frame/buffer_frames_x)
        segment_no_x = Int(round(frame/buffer_frames_x))
        MODT_aeif.save_rcdx_hdf5(string(rcd_name, "/all_x.h5"), segment_no_x, rcdx)
        MODT_aeif.initialize_rcdx!(rcdx)
        MODT_aeif.take_record_x!(0, system, rcdx, buffer_frames_x)
    end
end

end
