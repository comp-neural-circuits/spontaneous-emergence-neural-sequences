%Use fixed threshold 

early_net_dir = '../basic_network/2/';
net_dir = '../basic_network/3/';

if ~isfolder('strong_weights_turnover')
        mkdir('strong_weights_turnover')
end

resolution = 1000;

ckpt_num1 = 1001;     %total number of checkpoints

trials = 10;

thresh_EE = 2.5;
thresh_IE = 2.5;

V0 = -0.063;

start_point = 7100;   %starting point of new net
t_offset = 100;
early_t_offset = 100;
early_interval_w = 100;  %time interval to save weights for the previous net
       
all_ckpt = (ckpt_num1 - 4*resolution):resolution:ckpt_num1;
        
if all_ckpt(1) < 1
        plot_ckpt = (early_t_offset+1-start_point):early_interval_w:ckpt_num1;
else
        plot_ckpt = 1:ckpt_num1;
end

EE_ratio_strong = cell(length(all_ckpt), trials);
IE_ratio_strong = cell(length(all_ckpt), trials);

index_ee_strong = cell(length(all_ckpt), trials);
index_ie_strong = cell(length(all_ckpt), trials);

effw_EE_all = cell(length(all_ckpt), trials);
effw_IE_all = cell(length(all_ckpt), trials);

for trial = 1:trials
        folder = [net_dir, num2str(trial)];
        early_net_folder = [early_net_dir, num2str(trial)];
        for i = 1:length(all_ckpt)
                disp(['Resolution: ', num2str(resolution), 's, trial ' num2str(trial), ', checkpoint ', num2str(all_ckpt(i))])
                if all_ckpt(i) <= 1
                        EE_seg = h5read_sparsematrix([early_net_folder, '/all_w.h5'], ...
                                        [num2str((all_ckpt(i)+start_point-1)/early_interval_w), '/EE']);
                        IE_seg = h5read_sparsematrix([early_net_folder, '/all_w.h5'], ...
                                        [num2str((all_ckpt(i)+start_point-1)/early_interval_w), '/IE']);
                        ET = h5read([early_net_folder, '/all_w.h5'], ['/', num2str((all_ckpt(i)+start_point-1)/early_interval_w), '/ET']);
                        IT = h5read([early_net_folder, '/all_w.h5'], ['/', num2str((all_ckpt(i)+start_point-1)/early_interval_w), '/IT']);
                else
                        EE_seg = h5read_sparsematrix([folder, '/all_w.h5'], [num2str(all_ckpt(i) + t_offset - 1), '/EE']);
                        IE_seg = h5read_sparsematrix([folder, '/all_w.h5'], [num2str(all_ckpt(i) + t_offset -1), '/IE']);
                        ET = h5read([folder, '/all_w.h5'], ['/', num2str(all_ckpt(i) + t_offset - 1), '/ET']);
                        IT = h5read([folder, '/all_w.h5'], ['/', num2str(all_ckpt(i) + t_offset - 1), '/IT']);
                end
                
                th_e_post_ee = ET(EE_seg.rowval);
                th_e_post_ie = ET(IE_seg.rowval);
                
                effw_EE = EE_seg.nzval./(th_e_post_ee - V0)/1e3;
                effw_IE = IE_seg.nzval./(th_e_post_ie - V0)/1e3;
                
                effw_EE_all{i, trial} = effw_EE;
                effw_IE_all{i, trial} = effw_IE;

                index_ee_strong{i, trial} = find(effw_EE > thresh_EE);
                index_ie_strong{i, trial} = find(effw_IE > thresh_IE);
        
                EE_ratio_strong{i, trial} = zeros(length(plot_ckpt), length(index_ee_strong{i, trial}));

                IE_ratio_strong{i, trial} = zeros(length(plot_ckpt), length(index_ie_strong{i, trial}));
        
                for idx_ckpt = 1:length(plot_ckpt)
                        ckpt = plot_ckpt(idx_ckpt);
                        if ckpt <= 1
                                if ckpt == (1-start_point)
                                        EE_seg = h5read_sparsematrix([early_net_folder, '/init_condition_w.h5'], 'EE');
                                        IE_seg = h5read_sparsematrix([early_net_folder, '/init_condition_w.h5'], 'IE');
                                        ET = h5read([early_net_folder, '/init_condition_w.h5'], '/ET');
                                        IT = h5read([early_net_folder, '/init_condition_w.h5'], '/IT');
                                else
                                        EE_seg = h5read_sparsematrix([early_net_folder, '/all_w.h5'], ...
                                                        [num2str((ckpt+start_point-1)/early_interval_w), '/EE']);
                                        IE_seg = h5read_sparsematrix([early_net_folder, '/all_w.h5'], ...
                                                        [num2str((ckpt+start_point-1)/early_interval_w), '/IE']);
                                        ET = h5read([early_net_folder, '/all_w.h5'], ['/', num2str((ckpt+start_point-1)/early_interval_w), '/ET']);
                                        IT = h5read([early_net_folder, '/all_w.h5'], ['/', num2str((ckpt+start_point-1)/early_interval_w), '/IT']);
                                end
                        else
                                EE_seg = h5read_sparsematrix([folder, '/all_w.h5'], [num2str(ckpt + t_offset -1), '/EE']);
                                IE_seg = h5read_sparsematrix([folder, '/all_w.h5'], [num2str(ckpt + t_offset -1), '/IE']);
                                ET = h5read([folder, '/all_w.h5'], ['/', num2str(ckpt + t_offset -1), '/ET']);
                                IT = h5read([folder, '/all_w.h5'], ['/', num2str(ckpt + t_offset -1), '/IT']);
                        end
                        th_e_post_ee = ET(EE_seg.rowval);
                        th_e_post_ie = ET(IE_seg.rowval);
                        effw_EE = EE_seg.nzval./(th_e_post_ee - V0)/1e3;
                        effw_IE = IE_seg.nzval./(th_e_post_ie - V0)/1e3;
                        EE_ratio_strong{i, trial}(idx_ckpt, :) = effw_EE(index_ee_strong{i, trial})';
                        IE_ratio_strong{i, trial}(idx_ckpt, :) = effw_IE(index_ie_strong{i, trial})';
                end
        end
end

save('strong_weights_turnover/strong_weights_all_resolution.mat', 'EE_ratio_strong', 'IE_ratio_strong') 

init_net_dir = '../basic_network/1/';

init_EE_all = cell(1, trials);
init_IE_all = cell(1, trials);

for trial = 1:trials
        EE_seg = h5read_sparsematrix([init_net_dir, num2str(trial), '/init_condition_w.h5'], 'EE');
        IE_seg = h5read_sparsematrix([init_net_dir, num2str(trial), '/init_condition_w.h5'], 'IE');
        ET = h5read([init_net_dir, num2str(trial), '/init_condition_w.h5'], '/ET');
        effw_EE = EE_seg.nzval./(ET(EE_seg.rowval) - V0)/1e3;
        effw_IE = IE_seg.nzval./(ET(IE_seg.rowval) - V0)/1e3;
        init_EE_all{trial} = effw_EE;
        init_IE_all{trial} = effw_IE;
end
        
save('strong_weights_turnover/effw_distributions.mat', 'effw_EE_all', 'effw_IE_all', 'init_EE_all', 'init_IE_all')
        
        
        