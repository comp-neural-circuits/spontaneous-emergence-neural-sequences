using Statistics
using SparseArrays
using HDF5
include("run_process.jl")

@with_kw mutable struct Record_w
    frames::Int = Int(1e4)
    Ne::Int = 400
    Ni::Int = 80
    EE::SparseMatrixCSC{Float, Int} = zeros(Ne, Ne)
    IE::SparseMatrixCSC{Float, Int} = zeros(Ne, Ni)
    ET::Array{Float} = zeros(Ne)
    IT::Array{Float} = zeros(Ni)
end

function take_record_w!(frame::Int, system::the_whole_system, rcdw::Record_w, buffer_frames_w::Int)
    if mod(frame, buffer_frames_w) == 0
        rcdw.EE = sparse(system.EE.w)
        rcdw.IE = sparse(system.IE.w)
        rcdw.ET = system.E.T
        rcdw.IT = system.I.T
    end
end

function initialize_rcdw!(rcdw::Record_w)
    rcdw.EE = zeros(rcdw.Ne, rcdw.Ne)
    rcdw.IE = zeros(rcdw.Ne, rcdw.Ni)
    rcdw.ET = zeros(Ne)
    rcdw.IT = zeros(Ni)
end

function save_rcdw_hdf5(filename::String, segment_no::Int, rcdw::Record_w)
    h5write(filename, string(segment_no, "/frames"), rcdw.frames)
    h5write(filename, string(segment_no, "/Ne"), rcdw.Ne)
    h5write(filename, string(segment_no, "/Ni"), rcdw.Ni)
    h5write_sparsematrix(filename, string(segment_no, "/EE"), rcdw.EE)
    h5write_sparsematrix(filename, string(segment_no, "/IE"), rcdw.IE)
    h5write(filename, string(segment_no, "/ET"), rcdw.ET)
    h5write(filename, string(segment_no, "/IT"), rcdw.IT)
end

function h5write_sparsematrix(filename, varname, SM)
    h5write(filename, string(varname, "/m"), SM.m)
    h5write(filename, string(varname, "/n"), SM.n)
    h5write(filename, string(varname, "/colptr"), SM.colptr)
    h5write(filename, string(varname, "/rowval"), SM.rowval)
    h5write(filename, string(varname, "/nzval"), SM.nzval)
end

function record_w_init_cond(filename, system)
    h5write_sparsematrix(filename, "EE", sparse(system.EE.w))
    h5write_sparsematrix(filename, "IE", sparse(system.IE.w))
    h5write_sparsematrix(filename, "EI", sparse(system.EI.w))
    h5write_sparsematrix(filename, "II", sparse(system.II.w))
    h5write(filename, "ET", system.E.T)
    h5write(filename, "IT", system.I.T)
end
