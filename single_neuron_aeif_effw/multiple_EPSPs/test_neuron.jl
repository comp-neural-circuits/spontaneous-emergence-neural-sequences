using Main.MOD_single_neuron_aeif
using Statistics
using LinearAlgebra
using MAT
using Random
using SparseArrays
using HDF5

include("parameters.jl")

println("Setting up neurons and networks.")

T = 2
dt = 1e-4
frames = Int(round(T/dt))
interval_x = 1

E = MOD_single_neuron_aeif.NeuronArray(N = Ne)
E.param.μθ = μθe
E.param.σθ = σθe
E.param.V_reset = V_reset_e

system = MOD_single_neuron_aeif.the_whole_system(dt = dt, E = E)

rcd = MOD_single_neuron_aeif.Record(frames = frames, interval_x = interval_x, Ne = Ne)

MOD_single_neuron_aeif.initialize_test_system!(system)
MOD_single_neuron_aeif.take_record!(0, system, rcd)

for frame = 1:frames
    MOD_single_neuron_aeif.run_process!(system, rcd, frame)
end

include("convert2mat.jl")
file_mat = matopen(string("single_neuron_effw=", effw, ".mat"), "w")
write(file_mat, "rcd", rcd1)
write(file_mat, "ET", E.T)
close(file_mat)
