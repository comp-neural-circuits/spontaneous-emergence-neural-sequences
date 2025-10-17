using Main.MODT_aeif
using Statistics
using LinearAlgebra
using MAT
using Random
using SparseArrays
using HDF5

include("parameters.jl")
include("retrieve_sp.jl")

println("Setting up neurons and networks.")

T = 0.4
dt = 1e-4
frames = Int(round(T/dt))
interval_x = 1#Int(round(1/dt))
interval_w = Int(round(T/dt))

ckpt = 80

network_path = string("../../basic_network/1/1/")

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
EI.param.exci_inh = 1
EI.param.plastic = false
EI.param.delay = Int(round(d_ei/dt))

II = MODT_aeif.Network(pre = I, post = I)
II.w = II_w
II.sparse = sign.(II.w)
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
EE.param.exci_inh = 1
EE.param.plastic = false
EE.param.plasticity_rule = "eSTDP"
EE.param.delay = Int(round(d_ee/dt))
E.inject_point = ones(Ne)

IE = MODT_aeif.Network(pre = I, post = E)
if ckpt == 0
    IE.w = retrieve_sp(h5read(string(network_path, "init_condition_w.h5"), "IE"))
else
    IE.w = retrieve_sp(h5read(string(network_path, "all_w.h5"), string(ckpt, "/IE")))
end
IE.sparse = sign.(IE.w)
IE.param.exci_inh = -1
IE.param.plastic = false
IE.param.plasticity_rule = "iSTDP"
IE.param.delay = Int(round(d_ie/dt))

println("Checkpoint ", ckpt, ", allocating memory for records")

system = MODT_aeif.the_whole_system(dt = dt, E = E, I = I,
    EE = EE, EI = EI, IE = IE, II = II)

rcd = MODT_aeif.Record(frames = frames, interval_x = interval_x, interval_w =
    interval_w, Ne = Ne, Ni = Ni)

MODT_aeif.initialize_test_system!(system)
MODT_aeif.take_record!(0, system, rcd)

println("Checkpoint ", ckpt, ", start running simulation")

for frame = 1:frames
    if ceil(10*frame*dt) == floor(10*frame*dt)
        println("Checkpoint ", ckpt, ", ", frame*dt, " seconds")
    end
    MODT_aeif.run_process!(system, rcd, frame)
end

include("convert2mat.jl")
file_mat = matopen(string("testing_control.mat"), "w")
write(file_mat, "rcd", rcd1)
write(file_mat, "ET", E.T)
write(file_mat, "IT", I.T)
close(file_mat)
