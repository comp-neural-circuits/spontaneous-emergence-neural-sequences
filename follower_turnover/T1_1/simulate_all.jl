all_dir = readdir()

for j = 1:length(all_dir)
    if isdigit(all_dir[j][1]) && ~isdir(string(all_dir[j], "/sequence"))
        cd(all_dir[j])
        println(string("Runing test ", j))
        #include("start_simulation.jl")
        include(joinpath(pwd(), "start_simulation.jl"))
        cd("..")
    end
end
