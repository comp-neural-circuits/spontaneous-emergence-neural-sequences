using Statistics
using SparseArrays
include("run_process.jl")

@with_kw mutable struct Record
    frames::Int = Int(1e4)
    interval_x::Int = 10
    interval_w::Int = 100
    Ne::Int = 400
    Ni::Int = 80
    xe::SparseMatrixCSC{Float, Int} = zeros(Ne, Integer(round(frames/interval_x))+1)
    xi::SparseMatrixCSC{Float, Int} = zeros(Ni, Integer(round(frames/interval_x))+1)
    Ie::SparseMatrixCSC{Float, Int} = zeros(Ne, Integer(round(frames/interval_x))+1)
    Ii::SparseMatrixCSC{Float, Int} = zeros(Ni, Integer(round(frames/interval_x))+1)
    Ve::Array{Float}{2} = zeros(Ne, Int(round(frames/interval_x))+1)
    Vi::Array{Float}{2} = zeros(Ni, Int(round(frames/interval_x))+1)
    gee::Array{Float}{2} = zeros(Ne, Int(round(frames/interval_x))+1)
    gie::Array{Float}{2} = zeros(Ne, Int(round(frames/interval_x))+1)
    gei::Array{Float}{2} = zeros(Ni, Int(round(frames/interval_x))+1)
    gii::Array{Float}{2} = zeros(Ni, Int(round(frames/interval_x))+1)
    EE::Array{Any} = Array{Any}(undef, Int(round(frames/interval_w))+1)
    IE::Array{Any} = Array{Any}(undef, Int(round(frames/interval_w))+1)
end

function take_record!(frame::Int, system::the_whole_system, rcd::Record)
    if ceil(frame/rcd.interval_x) == floor(frame/rcd.interval_x)
        index_rcdV = Int(round(frame/rcd.interval_x) + 1)
        rcd.xe[:, index_rcdV] = sparse(system.E.x)
        rcd.xi[:, index_rcdV] = sparse(system.I.x)
        rcd.Ie[:, index_rcdV] = sparse(system.E.I)
        rcd.Ii[:, index_rcdV] = sparse(system.I.I)
        rcd.Ve[:, index_rcdV] = system.E.V
        rcd.Vi[:, index_rcdV] = system.I.V
        rcd.gee[:, index_rcdV] = system.E.ge
        rcd.gie[:, index_rcdV] = system.E.gi
        rcd.gei[:, index_rcdV] = system.I.ge
        rcd.gii[:, index_rcdV] = system.I.gi
    end
    if ceil(frame/rcd.interval_w) == floor(frame/rcd.interval_w)
        index_rcdw = Int(round((frame/rcd.interval_w)) + 1)
        rcd.EE[index_rcdw] = sparse(system.EE.w)
        rcd.IE[index_rcdw] = sparse(system.IE.w)
    end
end
