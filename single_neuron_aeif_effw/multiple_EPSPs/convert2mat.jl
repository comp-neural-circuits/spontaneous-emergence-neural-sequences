using SparseArrays
using Revise
using Parameters

@with_kw mutable struct Record_MAT
    frames::Int = Int(1e4)
    interval_x::Int = 10
    Ne::Int = 400
    Ni::Int = 0
    xe::SparseMatrixCSC{Float64, Int} = zeros(Ne, Integer(round(frames/interval_x))+1)
    Ve::Array{Float32}{2} = zeros(Ne, Int(round(frames/interval_x))+1)
    gee::Array{Float32}{2} = zeros(Ne, Int(round(frames/interval_x))+1)
    #=
    Ie::SparseMatrixCSC{Float64, Int} = zeros(Ne, Integer(round(frames/interval_x))+1)
    we::Array{Float}{2} = zeros(Ne, Int(round(frames/interval_x))+1)
    gie::Array{Float32}{2} = zeros(Ne, Int(round(frames/interval_x))+1)
    =#
end

#rcd1 = Record_MAT(frames = rcd.frames, interval_x = rcd.interval_x,
#Ne = rcd.Ne, Ve = rcd.Ve, gee = rcd.gee, gie = rcd.gie, we = rcd.we)

#rcd1 = Record_MAT(frames = rcd.frames, interval_x = rcd.interval_x, Ne = rcd.Ne)

rcd1 = Record_MAT(frames = rcd.frames, interval_x = rcd.interval_x,
Ne = rcd.Ne, Ve = rcd.Ve, gee = rcd.gee)

#rcd1.Ie = SparseMatrixCSC{Float64, Int}(rcd.Ie)
rcd1.xe = SparseMatrixCSC{Float64, Int}(rcd.xe)
