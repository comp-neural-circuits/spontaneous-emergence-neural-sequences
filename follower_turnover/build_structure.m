trials = 10;
source_neurons = cell(trials, 1);

source_neurons{1} = [1, 3, 5, 6, 9];
source_neurons{2} = [1, 3, 4, 5, 6];
source_neurons{3} = [1, 2, 3, 5, 6];
source_neurons{4} = [1, 2, 5, 6, 7];
source_neurons{5} = [1, 4, 6, 8, 9];
source_neurons{6} = [2, 3, 4, 5, 6];
source_neurons{7} = [1, 4, 5, 6, 7];
source_neurons{8} = [1, 3, 4, 5, 6];
source_neurons{9} = [1, 2, 4, 5, 6];
source_neurons{10} = [1, 2, 5, 6, 7];

template_dir = 'T1_1';

for i = 1:trials
        for j = 1:length(source_neurons{i})
                if i ~= 1 || j ~= 1 
                        new_dir = ['T', num2str(i), '_', num2str(source_neurons{i}(j))];
                        if ~isfolder(new_dir)
                                mkdir(new_dir)
                                mkdir([new_dir, '/template'])
                        end
                        
                        copyfile([template_dir, '/template'], [new_dir, '/template']);
                        copyfile([template_dir, '/build_structure.m'], new_dir);
                        copyfile([template_dir, '/delay_and_jitter_all.m'], new_dir);
                        copyfile([template_dir, '/find_seq_all.m'], new_dir);
                        copyfile([template_dir, '/replace_file.m'], new_dir);
                        copyfile([template_dir, '/simulate_all.jl'], new_dir);
                        
                        cd(new_dir);
                        build_structure;
                        cd('..');
                end
        end

end

