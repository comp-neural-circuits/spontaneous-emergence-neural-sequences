using SparseArrays
using Revise
using Parameters

@with_kw mutable struct Record_MAT
    frames::Int = Int(1e4)
    interval_x::Int = 10
    interval_w::Int = 100
    Ne::Int = 400
    Ni::Int = 80
    xe::SparseMatrixCSC{Float64, Int} = zeros(Ne, Integer(round(frames/interval_x))+1)
    xi::SparseMatrixCSC{Float64, Int} = zeros(Ni, Integer(round(frames/interval_x))+1)
    Ie::SparseMatrixCSC{Float64, Int} = zeros(Ne, Integer(round(frames/interval_x))+1)
    Ii::SparseMatrixCSC{Float64, Int} = zeros(Ni, Integer(round(frames/interval_x))+1)
    Ve::Array{Float32}{2} = zeros(Ne, Int(round(frames/interval_x))+1)
    Vi::Array{Float32}{2} = zeros(Ni, Int(round(frames/interval_x))+1)
    gee::Array{Float32}{2} = zeros(Ne, Int(round(frames/interval_x))+1)
    gie::Array{Float32}{2} = zeros(Ne, Int(round(frames/interval_x))+1)
    gei::Array{Float32}{2} = zeros(Ni, Int(round(frames/interval_x))+1)
    gii::Array{Float32}{2} = zeros(Ni, Int(round(frames/interval_x))+1)
    EE::Array{Any} = Array{Any}(undef, Int(round(frames/interval_w))+1)
    IE::Array{Any} = Array{Any}(undef, Int(round(frames/interval_w))+1)
end

rcd1 = Record_MAT(frames = rcd.frames, interval_x = rcd.interval_x,
interval_w = rcd.interval_w, Ne = rcd.Ne, Ni = rcd.Ni, Ve = rcd.Ve,
Vi = rcd.Vi, gee = rcd.gee, gie = rcd.gie, gei = rcd.gei, gii = rcd.gii)


rcd1.Ie = SparseMatrixCSC{Float64, Int}(rcd.Ie)
rcd1.Ii = SparseMatrixCSC{Float64, Int}(rcd.Ii)
rcd1.xe = SparseMatrixCSC{Float64, Int}(rcd.xe)
rcd1.xi = SparseMatrixCSC{Float64, Int}(rcd.xi)
for i = 1:length(rcd.EE)
    rcd1.EE[i] = SparseMatrixCSC{Float64, Int}(rcd.EE[i])
end
for i = 1:length(rcd.IE)
    rcd1.IE[i] = SparseMatrixCSC{Float64, Int}(rcd.IE[i])
end
