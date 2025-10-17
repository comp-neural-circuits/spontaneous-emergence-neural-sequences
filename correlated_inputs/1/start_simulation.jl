using Distributed

num_sim = 10

addprocs(num_sim)

try
    simulations = @distributed for sim_id = 1:num_sim
        include("MOD_aeif.jl")
        include("train_net.jl")
        train_net(sim_id)
    end
    fetch(simulations)
catch
    simulations = @distributed for sim_id = 1:num_sim
        include("MOD_aeif.jl")
        include("train_net.jl")
        train_net(sim_id)
    end
    fetch(simulations)
end
