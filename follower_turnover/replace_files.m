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

ckpt = [80:100, 110:10:160];

template_dir = 'T1_1';

for i = 1:trials
        for j = 1:length(source_neurons{i})
                new_dir = ['T', num2str(i), '_', num2str(source_neurons{i}(j))];
                        
                
                if i~=1 || j~=1
                %copyfile([template_dir, '/build_structure.m'], new_dir);
                copyfile([template_dir, '/delay_and_jitter_all.m'], new_dir);
                copyfile([template_dir, '/find_seq_all.m'], new_dir);
                %copyfile([template_dir, '/simulate_all.jl'], new_dir);
                end
                %}
                %{
                file_name = 'test_net.jl';
                if i~=1 || j~=1
                        copyfile([template_dir, '/template/', file_name], [new_dir, '/template/']);
                end
                for k = 1:length(ckpt)
                        copyfile([new_dir, '/template/', file_name], ...
                                [new_dir, '/', num2str(i), '_', num2str(source_neurons{i}(j)), '_', num2str(ckpt(k))]);
                end
                %}
        end

end