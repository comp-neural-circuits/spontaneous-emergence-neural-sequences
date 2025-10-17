using Statistics
using SparseArrays
include("run_process.jl")

@with_kw mutable struct Record
    frames::Int = Int(1e4)
    interval_x::Int = 10
    Ne::Int = 400
    xe::SparseMatrixCSC{Float, Int} = zeros(Ne, Integer(round(frames/interval_x))+1)
    Ve::Array{Float}{2} = zeros(Ne, Int(round(frames/interval_x))+1)
    gee::Array{Float}{2} = zeros(Ne, Int(round(frames/interval_x))+1)
    #=
    Ie::SparseMatrixCSC{Float, Int} = zeros(Ne, Integer(round(frames/interval_x))+1)
    we::Array{Float}{2} = zeros(Ne, Int(round(frames/interval_x))+1)
    gie::Array{Float}{2} = zeros(Ne, Int(round(frames/interval_x))+1)
    =#
end

function take_record!(frame::Int, system::the_whole_system, rcd::Record)
    if ceil(frame/rcd.interval_x) == floor(frame/rcd.interval_x)
        index_rcdV = Int(round(frame/rcd.interval_x) + 1)
        rcd.xe[:, index_rcdV] = sparse(system.E.x)
        rcd.Ve[:, index_rcdV] = system.E.V
        rcd.gee[:, index_rcdV] = system.E.ge
        #=
        rcd.Ie[:, index_rcdV] = sparse(system.E.I)
        rcd.we[:, index_rcdV] = system.E.w
        rcd.gie[:, index_rcdV] = system.E.gi
        =#
    end
end
