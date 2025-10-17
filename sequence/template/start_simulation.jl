using Distributed

num_sim = 10

addprocs(num_sim)

try
    simulations = @distributed for sim_id = 1:num_sim
        include("MODT_aeif.jl")
        include("test_net.jl")
        test_net(sim_id)
    end
    fetch(simulations)
catch
    simulations = @distributed for sim_id = 1:num_sim
        include("MODT_aeif.jl")
        include("test_net.jl")
        test_net(sim_id)
    end
    fetch(simulations)
end

all_procs = procs()
for i = 1:length(all_procs)
    if all_procs[i] != 1
        rmprocs(all_procs[i])
    end
end
