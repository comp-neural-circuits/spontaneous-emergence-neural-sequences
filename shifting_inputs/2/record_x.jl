using Statistics
using SparseArrays
using HDF5
include("run_process.jl")

@with_kw mutable struct Record_x
    frames::Int = Int(1e4)
    interval_x::Int = 10
    Ne::Int = 400
    Ni::Int = 80
    xe::SparseMatrixCSC{Float, Int} = zeros(Ne, Integer(round(frames/interval_x))+1)
    xi::SparseMatrixCSC{Float, Int} = zeros(Ni, Integer(round(frames/interval_x))+1)
end

function take_record_x!(frame::Int, system::the_whole_system, rcdx::Record_x, buffer_frames_x::Int)
    pseudo_frame = mod(frame, buffer_frames_x)
    if frame > 0 && pseudo_frame == 0
        pseudo_frame = buffer_frames_x
    end
    if ceil(pseudo_frame/rcdx.interval_x) == floor(pseudo_frame/rcdx.interval_x)
        index_rcdx = Int(round(pseudo_frame/rcdx.interval_x) + 1)
        rcdx.xe[:, index_rcdx] = sparse(system.E.x)
        rcdx.xi[:, index_rcdx] = sparse(system.I.x)
    end
end

function initialize_rcdx!(rcdx::Record_x)
    rcdx.xe = zeros(rcdx.Ne, Integer(round(rcdx.frames/rcdx.interval_x))+1)
    rcdx.xi = zeros(rcdx.Ni, Integer(round(rcdx.frames/rcdx.interval_x))+1)
end

function save_rcdx_hdf5(filename::String, segment_no::Int, rcdx::Record_x)
    h5write(filename, string(segment_no, "/frames"), rcdx.frames)
    h5write(filename, string(segment_no, "/interval_x"), rcdx.interval_x)
    h5write(filename, string(segment_no, "/Ne"), rcdx.Ne)
    h5write(filename, string(segment_no, "/Ni"), rcdx.Ni)
    h5write_sparsematrix(filename, string(segment_no, "/xe"), rcdx.xe)
    h5write_sparsematrix(filename, string(segment_no, "/xi"), rcdx.xi)
end

function h5write_sparsematrix(filename, varname, SM)
    h5write(filename, string(varname, "/m"), SM.m)
    h5write(filename, string(varname, "/n"), SM.n)
    h5write(filename, string(varname, "/colptr"), SM.colptr)
    h5write(filename, string(varname, "/rowval"), SM.rowval)
    h5write(filename, string(varname, "/nzval"), SM.nzval)
end
