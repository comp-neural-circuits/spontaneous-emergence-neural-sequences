using Statistics
using LinearAlgebra
using Random
include("parameters.jl")

@with_kw mutable struct NeuronParameter
    μθ = -50e-3
    σθ = 1e-3
    V_reset = -60e-3
end

@with_kw mutable struct NeuronArray
    param::NeuronParameter = NeuronParameter()
    N::Int = 100    #number of neurons
    hub_neurons::Array{Int} = zeros(N)
    I::Array{Float}{1} = zeros(N)   #external inputs
    I_net_e::Array{Float}{1} = zeros(N)   #inputs from other excitory neurons
    I_net_i::Array{Float}{1} = zeros(N)   #inputs from other inhibitory neurons
    V::Array{Float}{1} = zeros(N)   #membrane potential
    w::Array{Float}{1} = zeros(N)   #adaptation
    I_spike::Array{Float}{1} = zeros(N)
    x::Array{Bool}{1} = Array{Bool}(zeros(N))   #spikes
    T::Array{Float}{1} = zeros(N)   #firing thresholds
    ge::Array{Float}{1} = rand(N)*0.01   #excitatory channels
    gi::Array{Float}{1} = rand(N)*0.01   #inhibitory channels
    past15fires::Array{Float}{2} = zeros(N, 15) #record spikes in past 15 frames
end

function update_neurons!(G::NeuronArray, dt::Float, t::Float)
    if abs(t/wn_dt - round(t/wn_dt)) < 1e-3
        #G.I = G.rand_inputs[:, Int(floor(t/wn_dt))+1]
        G.I = broadcast(+, μ_wn, σ_wn*randn(G.N))
    end
    G.w = G.w + dt*(a_w*broadcast(-, G.V, El) - G.w)/τw
    G.I_spike = gL * ΔT * exp.((G.V - G.T)/ΔT)
    G.V = G.V + dt*(gL*broadcast(-, El, G.V) + G.ge .* broadcast(-, Ee, G.V) +
    G.gi .* broadcast(-, Ei, G.V) + G.I + G.I_spike - G.w)/Cm
#=
    if abs((t-tstim)/dt) < 1e-3
        G.V = G.V + broadcast(-, Vpeak+1e-4, G.V) .* G.hub_neurons
    end
=#
    G.x = (G.V .> Vpeak)
    G.V[G.x] = G.param.V_reset * ones(Int(sum(G.x)))
    G.w[G.x] = broadcast(+, b_w, G.w[G.x])
    G.past15fires = [G.past15fires[:, 2:15]'; G.x']'
    G.ge = G.ge + G.I_net_e
    G.ge = G.ge - dt * G.ge/τe
    G.gi = G.gi + G.I_net_i
    G.gi = G.gi - dt * G.gi/τi
    if abs(t-tstim) < 1e-3*dt
        G.ge = G.ge + effw * broadcast(-, G.T, V0) * 1e3
    end
end

function initialize_neurons!(G::NeuronArray)
    G.I = broadcast(+, μ_wn, σ_wn*randn(G.N))
    G.I_net_e = zeros(G.N)
    G.I_net_i = zeros(G.N)
    G.V = broadcast(+, El+μ_wn/(a_w+gL), 0e-3*rand(G.N))
    G.w = a_w*broadcast(-, G.V, El)
    G.T = broadcast(+, G.param.σθ*randn(G.N), G.param.μθ)
    G.x = (G.V .> G.T)
    G.ge = zeros(G.N)
    G.gi = zeros(G.N)
end

function initialize_test_neurons!(G::NeuronArray)
    G.I = broadcast(+, μ_wn, σ_wn*randn(G.N))
    G.I_net_e = zeros(G.N)
    G.I_net_i = zeros(G.N)
    G.V = broadcast(+, El+μ_wn/(a_w+gL), 0e-3*rand(G.N))
    G.w = a_w*broadcast(-, G.V, El)
    G.x = (G.V .> G.T)
    G.ge = zeros(G.N)
    G.gi = zeros(G.N)
end
