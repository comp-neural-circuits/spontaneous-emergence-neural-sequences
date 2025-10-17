using Statistics
using LinearAlgebra
using Random
include("parameters.jl")

@with_kw mutable struct NeuronParameter
    μθ = -50e-3
    σθ = 1e-3
    V_reset = -60e-3
    rt = 0.5    #target rate
    ηIP = 0.05e-3   #update rate of IP
    p_corr = 1.0    #correlation, p_corr^2
    corr_t = 10.0
    gap_t = 5.0
    ξext = 5e-3
    τext = 40e-3    #poisson inputs
    time_offset = 8e3   #starting time
end

@with_kw mutable struct NeuronArray
    param::NeuronParameter = NeuronParameter()
    N::Int = 100    #number of neurons
    M::Int = 12      #number of correlated groups
    group_size::Int = 100   #number of neurons in each correlated groups
    groups::Array{Int}{2} = reshape(1:(M*group_size), M, group_size)'
    I::Array{Float}{1} = zeros(N)   #external inputs
    I_source::Float = 0   #source poisson inputs for correlation
    I_poisson::Array{Float}{1} = zeros(N)   #poisson inputs
    I_net_e::Array{Float}{1} = zeros(N)   #inputs from other excitory neurons
    I_net_i::Array{Float}{1} = zeros(N)   #inputs from other inhibitory neurons
    V::Array{Float}{1} = zeros(N)   #membrane potential
    w::Array{Float}{1} = zeros(N)   #adaptation
    I_spike::Array{Float}{1} = zeros(N)
    x::Array{Bool}{1} = Array{Bool}(zeros(N))   #spikes
    T::Array{Float64}{1} = zeros(N)   #firing thresholds
    ge::Array{Float}{1} = rand(N)*0.01   #excitatory channels
    gi::Array{Float}{1} = rand(N)*0.01   #inhibitory channels
    pastfires::Array{Float}{2} = zeros(N, 15) #record spikes in past 15 frames
end

function update_neurons!(G::NeuronArray, dt::Float, t::Float)
    if abs(t/wn_dt - round(t/wn_dt)) < 1e-3
        G.I = broadcast(+, μ_wn, σ_wn*randn(G.N))
    end
    if G.M > 0
        G.I_source = (rand() .< dt/G.param.τext)
        corr_period = Int(ceil((t + G.param.time_offset)/(G.param.corr_t + G.param.gap_t)))
        corr_idx = mod((corr_period-1), G.M)+1
        if t + G.param.time_offset - (corr_period-1)*(G.param.corr_t + G.param.gap_t) <= G.param.corr_t
            G.I_poisson[G.groups[:, corr_idx]] = G.param.ξext * Cm/dt *
            ((rand(G.group_size) .< (1-G.param.p_corr)*dt/G.param.τext) +
            G.I_source * (rand(G.group_size) .< G.param.p_corr))
            for i = 1:G.M
                if i != corr_idx
                    G.I_poisson[G.groups[:, i]] = G.param.ξext * Cm/dt * (rand(G.group_size) .< dt/G.param.τext)
                end
            end
        else
            G.I_poisson = G.param.ξext * Cm/dt * (rand(G.N) .< dt/G.param.τext)
        end
    end
    G.I_poisson[(G.M*G.group_size+1):end] = G.param.ξext * Cm/dt * (rand(G.N - G.M*G.group_size, 1) .< dt/G.param.τext)
    G.ge = G.ge + G.I_net_e
    G.gi = G.gi + G.I_net_i
    G.w = G.w + dt*(a_w*broadcast(-, G.V, El) - G.w)/τw
    G.I_spike = gL * ΔT * exp.((G.V - G.T)/ΔT)
    G.V = G.V + dt*(gL*broadcast(-, El, G.V) + G.ge .* broadcast(-, Ee, G.V) +
    G.gi .* broadcast(-, Ei, G.V) - G.w + G.I + G.I_spike + G.I_poisson)/Cm
    G.x = (G.V .> Vpeak)
    G.V[G.x] .= G.param.V_reset
    G.w[G.x] = broadcast(+, b_w, G.w[G.x])
    G.pastfires = [G.pastfires[:, 2:end]'; G.x']'
    G.ge = G.ge - dt * G.ge/τe
    G.gi = G.gi - dt * G.gi/τi
    G.T = G.T + Array{Float64}(G.param.ηIP * broadcast(-, G.x, G.param.rt*dt))
end

function initialize_neurons!(G::NeuronArray)
    G.I = broadcast(+, μ_wn, σ_wn*randn(G.N))
    G.I_poisson = zeros(G.N)
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
    G.I_poisson = zeros(G.N)
    G.I_net_e = zeros(G.N)
    G.I_net_i = zeros(G.N)
    G.V = broadcast(+, El+μ_wn/(a_w+gL), 0e-3*rand(G.N))
    G.w = a_w*broadcast(-, G.V, El)
    G.x = (G.V .> G.T)
    G.ge = zeros(G.N)
    G.gi = zeros(G.N)
end
