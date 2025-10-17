module MOD_aeif
using Revise
using Parameters
const Int = Int64
const Float = Float32

include("neurons.jl")
include("networks.jl")
include("record_x.jl")
include("record_w.jl")
include("run_process.jl")
include("parameters.jl")

end
