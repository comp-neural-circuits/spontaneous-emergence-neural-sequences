%decouple strong weights switching and hub neuron switching
%Use fixed ratio effw = 2.5
%Only E-to-E

net_dir = '../shifting_inputs/2/';

if ~isfolder('strong_weights_turnover')
        mkdir('strong_weights_turnover')
end

resolution = [100, 200, 500, 1000, 1500, 2000, 3000, 4000, 5000, 6000, 7000, 8000];

ckpt_num1 = 81;     %total number of checkpoints

thresh_EE = 2.5;

V0 = -0.063;
Ne = 1200;

trials = 10;

start_point = 8000;   %starting point of new net
t_offset = 100;
interval_w = 100;  %time interval to save weights for the new net

corr_range = 100;

n0_ee = cell(trials, length(resolution)); %remain
n1_ee = cell(trials, length(resolution)); %turnover

seg_num = round(Ne/corr_range);

parfor r = 1:length(resolution)
        
        all_ckpt = max(1, (ckpt_num1 - 10*round(resolution(r)/interval_w))): ...
                round(resolution(r)/interval_w):ckpt_num1;

        index_ee_strong = cell(length(all_ckpt), seg_num);

        for trial = 1:trials
                folder = [net_dir, num2str(trial)];
                for i = 1:length(all_ckpt)
                        disp(['Resolution: ', num2str(resolution(r)), 's, trial ' num2str(trial), ', checkpoint ', num2str(all_ckpt(i))])
                        
                        EE_seg = h5read_sparsematrix([folder, '/all_w.h5'], [num2str(all_ckpt(i)), '/EE']);
                        ET = h5read([folder, '/all_w.h5'], ['/', num2str(all_ckpt(i)), '/ET']);
                        
                        th_e_post_ee = ET(EE_seg.rowval);
                
                        effw_EE = EE_seg.nzval./(th_e_post_ee - V0)/1e3;
                        
                        EE_colval = zeros(length(EE_seg.nzval), 1);
                        for neuron = 1:EE_seg.n
                                index = (EE_seg.colptr(neuron) : (EE_seg.colptr(neuron+1) - 1));
                                EE_colval(index) = neuron;
                        end
                        
                        for corr_seg = 1:seg_num
                                index_ee_strong{i, corr_seg} = find(effw_EE > thresh_EE & EE_seg.rowval >= ...
                                        (corr_seg - 1)*corr_range + 1 & EE_seg.rowval <= corr_seg*corr_range & ...
                                        EE_colval >= (corr_seg - 1)*corr_range + 1 & EE_colval <= corr_seg*corr_range);
                        end
                       
                end

                n0_ee{trial, r} = zeros(length(all_ckpt) - 1, seg_num);     %no change
                n1_ee{trial, r} = zeros(length(all_ckpt) - 1, seg_num);     %turnover
                
                for i = 1:(length(all_ckpt) - 1)
                        for corr_seg = 1:seg_num
                                n0_ee{trial, r}(i, corr_seg) = min(length(index_ee_strong{i+1, corr_seg}), length(index_ee_strong{i, corr_seg}));
                                n1_ee{trial, r}(i, corr_seg) = length(index_ee_strong{i, corr_seg}) - n0_ee{trial, r}(i, corr_seg);
                        end
                end
                
        end
end

save('strong_weights_turnover/turnover_corr_diag_block_rate.mat', 'n0_ee', 'n1_ee', 'resolution')



 
