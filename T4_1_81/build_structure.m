groups = 19;
group_width = 0.5;
samples = 10;

V0 = -0.063;

EE_dir = '../../shifting_inputs/2/1/';
ckpt = 81;
EE_seg = h5read_sparsematrix([EE_dir, '/all_w.h5'], [num2str(ckpt), '/EE']);
EE = retrieve_sp(EE_seg);
ET = h5read([EE_dir, '/all_w.h5'], ['/', num2str(ckpt), '/ET']);
effw = max(EE./((ET - V0)*ones(1,size(EE, 2)))/1e3);

temp_dir = 'template';

for i = 1:groups
        mkdir(num2str(i));
        effw_start = (i-1)*group_width;
        effw_end = i*group_width;
        if i == groups
                neuron_idx = find(effw >= effw_start);
        else
                neuron_idx = find((effw >= effw_start) & (effw < effw_end));
        end
        if length(neuron_idx) <= samples
                neuron_idx_select = neuron_idx;
        else
                neuron_idx_select = neuron_idx(randperm(length(neuron_idx), samples));
        end
        cd(num2str(i));
        for j = 1:length(neuron_idx_select)
                mkdir(num2str(neuron_idx_select(j)));
                copyfile(['../', temp_dir], num2str(neuron_idx_select(j)));
        end
        cd('..');
end


        