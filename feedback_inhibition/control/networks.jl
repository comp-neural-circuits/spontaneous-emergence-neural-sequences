using Statistics
using SparseArrays
include("parameters.jl")

@with_kw mutable struct NetworkParameter
    plastic::Bool = true
    exci_inh::Int = 1
    plasticity_rule = "eSTDP"
    delay::Int = 15  #delayed frames
    lowlimit::Float = 1e-20
    A1::Float = 0
    A2::Float = 0   #amplitudes of synaptic modification
    τ1::Float = 1
    τ2::Float = 1   #time scales of A1 and A2
    LTDα::Float = 0
end

@with_kw mutable struct Network
    param::NetworkParameter = NetworkParameter()
    pre::NeuronArray
    post::NeuronArray
    w::SparseMatrixCSC{Float,Int64} = sprand(Float, post.N, pre.N, 1.0)
    wt::SparseMatrixCSC{Float,Int64} = SparseMatrixCSC(w')
    sparse::SparseMatrixCSC{Float,Int64} = sprand(Bool, post.N, pre.N, 1.0)
    outputs::Array{Float}{1} = zeros(post.N)
    Pa_pre::Array{Float}{1} = zeros(pre.N)
    Pa_post::Array{Float}{1} = zeros(post.N)
    sum_post::Array{Float}{1} = zeros(post.N)   #for Normalization
    sum_post_fixed::Array{Float}{1} = zeros(post.N)   #for Normalization
    tight_postSN::Bool = false
    sum_pre::Array{Float}{1} = zeros(pre.N)   #for Normalization
    sum_pre_fixed::Array{Float}{1} = zeros(pre.N)   #for Normalization
end

function update_output!(GG::Network, dt::Float)
    new_outputs = zeros(GG.post.N)
    for pre_idx = 1:GG.w.n
        if GG.pre.pastfires[pre_idx, end+1-GG.param.delay] == 1
            for syn_idx = GG.w.colptr[pre_idx]:(GG.w.colptr[pre_idx+1]-1)
                post_idx = GG.sparse.rowval[syn_idx]
                new_outputs[post_idx] = new_outputs[post_idx] + GG.w.nzval[syn_idx]
            end
        end
    end
    if GG.param.exci_inh == 1
        GG.post.I_net_e = GG.post.I_net_e - GG.outputs + new_outputs
    elseif GG.param.exci_inh == -1
        GG.post.I_net_i = GG.post.I_net_i - GG.outputs + new_outputs
    end
    GG.outputs = copy(new_outputs)
end

function update_network!(GG::Network, dt::Float, t::Float)
    if abs((t-tstim)/dt) < 1e-3
        GG.w .= 0
    end
    if GG.param.plastic
        GG.Pa_pre = GG.Pa_pre + GG.param.A1 * (GG.pre.x)
        GG.Pa_pre = GG.Pa_pre - dt*GG.Pa_pre/GG.param.τ1
        GG.Pa_post = GG.Pa_post + GG.param.A2 * (GG.post.x)
        GG.Pa_post = GG.Pa_post - dt*GG.Pa_post/GG.param.τ2
        for pre_idx = 1:GG.w.n
            if GG.pre.x[pre_idx] == 1
                for syn_idx = GG.w.colptr[pre_idx]:(GG.w.colptr[pre_idx+1]-1)
                    post_idx = GG.sparse.rowval[syn_idx]
                    w_old = GG.w[post_idx, pre_idx]
                    w_new = max(GG.param.lowlimit, w_old + GG.Pa_post[post_idx] - GG.param.LTDα)
                    GG.w[post_idx, pre_idx] = w_new
                    GG.wt[pre_idx, post_idx] = w_new
                    GG.sum_post[post_idx] = GG.sum_post[post_idx] + w_new - w_old
                    GG.sum_pre[pre_idx] = GG.sum_pre[pre_idx] + w_new - w_old
                end
            end
        end
        for post_idx = 1:GG.w.m
            if GG.post.x[post_idx] == 1
                for syn_idx = GG.wt.colptr[post_idx]:(GG.wt.colptr[post_idx+1]-1)
                    pre_idx = GG.wt.rowval[syn_idx]
                    w_old = GG.w[post_idx, pre_idx]
                    w_new = max(GG.param.lowlimit, w_old + GG.Pa_pre[pre_idx])
                    GG.w[post_idx, pre_idx] = w_new
                    GG.wt[pre_idx, post_idx] = w_new
                    GG.sum_post[post_idx] = GG.sum_post[post_idx] + w_new - w_old
                    GG.sum_pre[pre_idx] = GG.sum_pre[pre_idx] + w_new - w_old
                end
            end
        end
        preSyn_Norm!(GG)
        postSyn_Norm!(GG)
    end
    update_output!(GG, dt)
end

function initialize_network!(GG::Network)
    GG.outputs = zeros(GG.post.N)
end

function postSyn_Norm!(GG::Network)
    ratio = GG.sum_post_fixed ./ GG.sum_post
    if maximum(ratio) > 1+min_SN_correction || minimum(ratio) < 1-min_SN_correction
        for i = 1:length(GG.wt.nzval)
            GG.w.nzval[i] = GG.w.nzval[i]*ratio[GG.w.rowval[i]]
        end
        GG.wt = SparseMatrixCSC(GG.w')
        GG.sum_post = copy(GG.sum_post_fixed)
        for i = 1:GG.pre.N
            GG.sum_pre[i] = sum(GG.w.nzval[GG.w.colptr[i]:(GG.w.colptr[i+1]-1)])
        end
    end
end

function preSyn_Norm!(GG::Network)
    ratio = GG.sum_pre_fixed ./ GG.sum_pre
    if maximum(ratio) > 1+min_SN_correction || minimum(ratio) < 1-min_SN_correction
        for pre_idx = 1:GG.w.n
            for syn_idx = GG.w.colptr[pre_idx]:(GG.w.colptr[pre_idx+1]-1)
                GG.w.nzval[syn_idx] = GG.w.nzval[syn_idx]*ratio[pre_idx]
            end
        end
        GG.wt = SparseMatrixCSC(GG.w')
        GG.sum_pre = copy(GG.sum_pre_fixed)
        for i = 1:GG.post.N
            GG.sum_post[i] = sum(GG.wt.nzval[GG.wt.colptr[i]:(GG.wt.colptr[i+1]-1)])
        end
    end
end
