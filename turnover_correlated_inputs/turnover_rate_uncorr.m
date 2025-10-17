%decouple strong weights switching and hub neuron switching
%Use fixed ratio effw = 2.5

net_dir = '../correlated_inputs/2/';

if ~isfolder('strong_weights_turnover')
        mkdir('strong_weights_turnover')
end

resolution = [100, 200, 500, 1000, 1500, 2000, 3000, 4000, 5000, 6000, 7000, 8000];

ckpt_num1 = 151;     %total number of checkpoints

thresh_EE = 2.5;
thresh_IE = 2.5;

V0 = -0.063;
Ne = 1200;
Ni = 240;

trials = 10;

start_point = 10000;   %starting point of new net
t_offset = 100;
interval_w = 100;  %time interval to save weights for the new net

corr_range = 100;

n0_ee = cell(trials, length(resolution));   %change
n1_ee = cell(trials, length(resolution));   %no change
n0_ie = cell(trials, length(resolution));
n1_ie = cell(trials, length(resolution));

parfor r = 1:length(resolution)
        
        all_ckpt = max(1, (ckpt_num1 - 10*round(resolution(r)/interval_w))): ...
                round(resolution(r)/interval_w):ckpt_num1;
        
        index_ee_strong = cell(length(all_ckpt), 1);
        index_ie_strong = cell(length(all_ckpt), 1);

        for trial = 1:trials
                folder = [net_dir, num2str(trial)];
                for i = 1:length(all_ckpt)
                        disp(['Resolution: ', num2str(resolution(r)), 's, trial ' num2str(trial), ', checkpoint ', num2str(all_ckpt(i))])
                        
                        EE_seg = h5read_sparsematrix([folder, '/all_w.h5'], [num2str(all_ckpt(i)), '/EE']);
                        IE_seg = h5read_sparsematrix([folder, '/all_w.h5'], [num2str(all_ckpt(i)), '/IE']);
                        ET = h5read([folder, '/all_w.h5'], ['/', num2str(all_ckpt(i)), '/ET']);
                        IT = h5read([folder, '/all_w.h5'], ['/', num2str(all_ckpt(i)), '/IT']);
                
                        th_e_post_ee = ET(EE_seg.rowval);
                        th_e_post_ie = ET(IE_seg.rowval);
                
                        effw_EE = EE_seg.nzval./(th_e_post_ee - V0)/1e3;
                        effw_IE = IE_seg.nzval./(th_e_post_ie - V0)/1e3;
                        
                        EE_colval = zeros(length(EE_seg.nzval), 1);
                        for neuron = 1:EE_seg.n
                                index = (EE_seg.colptr(neuron) : (EE_seg.colptr(neuron+1) - 1));
                                EE_colval(index) = neuron;
                        end
                        
                        IE_colval = zeros(length(IE_seg.nzval), 1);
                        for neuron = 1:IE_seg.n
                                index = (IE_seg.colptr(neuron) : (IE_seg.colptr(neuron+1) - 1));
                                IE_colval(index) = neuron;
                        end

                        index_ee_strong{i} = [];
                        for corr_seg = 1:round(EE_seg.m/corr_range)
                                index_ee_strong{i} = [index_ee_strong{i};
                                find(effw_EE > thresh_EE & (EE_seg.rowval <= (corr_seg - 1)*corr_range | ...
                                EE_seg.rowval > corr_seg*corr_range) & EE_colval >= (corr_seg - 1)*corr_range + 1 & ...
                                EE_colval <= corr_seg*corr_range)];
                        end
                        index_ie_strong{i} = find(effw_IE > thresh_IE);
                end

                n0_ee{trial, r} = zeros(length(all_ckpt) - 1, 1);     %no change
                n1_ee{trial, r} = zeros(length(all_ckpt) - 1, 1);     %change
                for i = 1:(length(all_ckpt) - 1)
                        n0_ee{trial, r}(i) = length(intersect(index_ee_strong{i}, index_ee_strong{i+1}));
                        n1_ee{trial, r}(i) = length(index_ee_strong{i}) - n0_ee{trial, r}(i);
                end
                
                n0_ie{trial, r} = zeros(length(all_ckpt) - 1, 1);     %no change
                n1_ie{trial, r} = zeros(length(all_ckpt) - 1, 1);     %change
                for i = 1:(length(all_ckpt) - 1)
                        n0_ie{trial, r}(i) = length(intersect(index_ie_strong{i}, index_ie_strong{i+1}));
                        n1_ie{trial, r}(i) = length(index_ie_strong{i}) - n0_ie{trial, r}(i);
                end
                
        end
end

ratio_n0_n1 = zeros(length(resolution), trials, 2, 2);

for r = 1:length(resolution)
        n_ee_all = sum(cell2mat(n0_ee(:, r)') + cell2mat(n1_ee(:, r)'), 1);
        n_ie_all = sum(cell2mat(n0_ie(:, r)') + cell2mat(n1_ie(:, r)'), 1);
        ratio_n0_n1(r, :, 1, 1) = sum(cell2mat(n0_ee(:, r)'), 1) ./ n_ee_all;
        ratio_n0_n1(r, :, 1, 2) = sum(cell2mat(n1_ee(:, r)'), 1) ./ n_ee_all;
        ratio_n0_n1(r, :, 2, 1) = sum(cell2mat(n0_ie(:, r)'), 1) ./ n_ie_all;
        ratio_n0_n1(r, :, 2, 2) = sum(cell2mat(n1_ie(:, r)'), 1) ./ n_ie_all;
end

save('strong_weights_turnover/turnover_uncorr_rate.mat', 'ratio_n0_n1', 'n0_ee', 'n1_ee', 'n0_ie', 'n1_ie', 'resolution')



 
