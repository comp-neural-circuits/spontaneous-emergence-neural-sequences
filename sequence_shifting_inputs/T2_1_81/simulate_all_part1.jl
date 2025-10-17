group_inds = 1:5
groups = length(group_inds)

for i = 1:groups
    cd(string(group_inds[i]))
    all_dir = readdir()
    for j = 1:length(all_dir)
        cd(all_dir[j])
        println(string("Runing group ", group_inds[i], ", neuron ", j))
        #include("start_simulation.jl")
        include(joinpath(pwd(), "start_simulation.jl"))
        cd("..")
    end
    cd("..")
end
