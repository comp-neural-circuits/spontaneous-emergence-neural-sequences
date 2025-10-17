module MOD_single_neuron_aeif
using Revise
using Parameters
const Int = Int64
const Float = Float32

include("neurons.jl")
include("record.jl")
include("run_process.jl")
include("parameters.jl")

end
