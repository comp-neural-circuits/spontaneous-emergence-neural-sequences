groups = 1;

neuron_idx{1} = [1, 68, 80, 87, 116, 117, 420, 448, 546, 591, 669, 680, 752, ...
        797, 1006, 1013, 1065, 1081, 1196];

temp_dir = 'template';

for i = 1:groups
        mkdir(num2str(i));
        
        cd(num2str(i));
        for j = 1:length(neuron_idx{i})
                mkdir(num2str(neuron_idx{i}(j)));
                copyfile(['../', temp_dir], num2str(neuron_idx{i}(j)));
        end
        cd('..');
end


        