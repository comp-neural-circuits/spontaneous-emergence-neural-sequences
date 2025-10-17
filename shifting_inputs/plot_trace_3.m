%6(5000), 7(3000), 8(8000)

dir_all = {'1', '4'};
duration_all = [80, 80]*1e2;
t_start = [0, 100];
reorder = [0, 1];

trials = 10;

w_interval = 100;

V0 = -0.063;

th_0 = 2.5;  %threshold for high weight-thresh ratio

Ne = 1200;

p_ee = 0.17/3;

corr_range = 100;

corr_seg = round(Ne/corr_range);

ckpt_num_all = round(sum(duration_all)/w_interval)+1;

if ~isfolder('trace_3')
        mkdir('trace_3')
end
%{
avg_EE_corr_diag = cell(length(dir_all), trials);
avg_EE_corr_subdiag = cell(length(dir_all), trials);
avg_EE_uncorr = cell(length(dir_all), trials);
std_EE_corr_diag = cell(length(dir_all), trials);
std_EE_corr_subdiag = cell(length(dir_all), trials);
std_EE_uncorr = cell(length(dir_all), trials);
low_EE_corr_diag = cell(length(dir_all), trials);
low_EE_corr_subdiag = cell(length(dir_all), trials);
high_EE_corr_diag = cell(length(dir_all), trials);
high_EE_corr_subdiag = cell(length(dir_all), trials);
low_EE_uncorr = cell(length(dir_all), trials);
high_EE_uncorr = cell(length(dir_all), trials);
avg_EE_all = cell(length(dir_all), trials);

avg_EE_corr_diag_reorder = cell(length(dir_all), trials);
avg_EE_corr_subdiag_reorder = cell(length(dir_all), trials);
avg_EE_uncorr_reorder = cell(length(dir_all), trials);
std_EE_corr_diag_reorder = cell(length(dir_all), trials);
std_EE_corr_subdiag_reorder = cell(length(dir_all), trials);
std_EE_uncorr_reorder = cell(length(dir_all), trials);
low_EE_corr_diag_reorder = cell(length(dir_all), trials);
low_EE_corr_subdiag_reorder = cell(length(dir_all), trials);
high_EE_corr_diag_reorder = cell(length(dir_all), trials);
high_EE_corr_subdiag_reorder = cell(length(dir_all), trials);
low_EE_uncorr_reorder = cell(length(dir_all), trials);
high_EE_uncorr_reorder = cell(length(dir_all), trials);
avg_EE_all_reorder = cell(length(dir_all), trials);

for dir = 1:length(dir_all)
        for trial = 1:trials
                net_dir = [dir_all{dir}, '/', num2str(trial), '/'];
                ckpt_num = round(duration_all(dir)/w_interval) + 1;
                ckpt_start = round(t_start(dir)/w_interval) + 1;

                avg_EE_corr_diag{dir, trial} = zeros(ckpt_num, 1);
                avg_EE_corr_subdiag{dir, trial} = zeros(ckpt_num, 1);
                avg_EE_uncorr{dir, trial} = zeros(ckpt_num, 1);
                std_EE_corr_diag{dir, trial} = zeros(ckpt_num, 1);
                std_EE_corr_subdiag{dir, trial} = zeros(ckpt_num, 1);
                std_EE_uncorr{dir, trial} = zeros(ckpt_num, 1);
                low_EE_corr_diag{dir, trial} = zeros(ckpt_num, 1);
                low_EE_corr_subdiag{dir, trial} = zeros(ckpt_num, 1);
                high_EE_corr_diag{dir, trial} = zeros(ckpt_num, 1);
                high_EE_corr_subdiag{dir, trial} = zeros(ckpt_num, 1);
                low_EE_uncorr{dir, trial} = zeros(ckpt_num, 1);
                high_EE_uncorr{dir, trial} = zeros(ckpt_num, 1);

                avg_EE_corr_diag_reorder{dir, trial} = zeros(ckpt_num, 1);
                avg_EE_corr_subdiag_reorder{dir, trial} = zeros(ckpt_num, 1);
                avg_EE_uncorr_reorder{dir, trial} = zeros(ckpt_num, 1);
                std_EE_corr_diag_reorder{dir, trial} = zeros(ckpt_num, 1);
                std_EE_corr_subdiag_reorder{dir, trial} = zeros(ckpt_num, 1);
                std_EE_uncorr_reorder{dir, trial} = zeros(ckpt_num, 1);
                low_EE_corr_diag_reorder{dir, trial} = zeros(ckpt_num, 1);
                low_EE_corr_subdiag_reorder{dir, trial} = zeros(ckpt_num, 1);
                high_EE_corr_diag_reorder{dir, trial} = zeros(ckpt_num, 1);
                high_EE_corr_subdiag_reorder{dir, trial} = zeros(ckpt_num, 1);
                low_EE_uncorr_reorder{dir, trial} = zeros(ckpt_num, 1);
                high_EE_uncorr_reorder{dir, trial} = zeros(ckpt_num, 1);

                avg_EE_all{dir, trial} = zeros(corr_seg, corr_seg, ckpt_num);
                
                for ckpt = 1:ckpt_num
                        if ckpt + ckpt_start - 1 == 1
                                EE_seg = h5read_sparsematrix([net_dir, '/init_condition_w.h5'], 'EE');
                                ET = h5read([net_dir, '/init_condition_w.h5'], '/ET');
                        else
                                EE_seg = h5read_sparsematrix([net_dir, 'all_w.h5'], [num2str(ckpt + ckpt_start - 2), '/EE']);
                                ET = h5read([net_dir, 'all_w.h5'], ['/', num2str(ckpt + ckpt_start - 2), '/ET']);
                        end
                        EE = retrieve_sp(EE_seg);
                        effw = EE./(ET*ones(1, Ne) - V0)/1e3;
                        avg_EE = zeros(corr_seg, corr_seg);
                        reorder_idx = reshape(reshape(1:Ne, round(Ne/corr_range), corr_range)', 1, Ne);
                        effw_reorder = effw(reorder_idx, reorder_idx);
                        avg_EE_reorder = zeros(corr_seg, corr_seg);
                        for i = 1:corr_seg
                                for j = 1:corr_seg
                                        avg_EE(i, j) = mean(mean(effw((1:corr_range)+(i-1)*corr_range, ...
                                                (1:corr_range)+(j-1)*corr_range) > th_0))/p_ee;
                                        avg_EE_reorder(i, j) = mean(mean(effw_reorder((1:corr_range)+(i-1)*corr_range, ...
                                                       (1:corr_range)+(j-1)*corr_range) > th_0))/p_ee;
                                end
                        end
                        avg_EE_all{dir, trial}(:, :, ckpt) = avg_EE;
                        avg_EE_corr_diag{dir, trial}(ckpt) = mean(diag(avg_EE));
                        std_EE_corr_diag{dir, trial}(ckpt) = std(diag(avg_EE), 1);
                        low_EE_corr_diag{dir, trial}(ckpt) = prctile(diag(avg_EE), 0);
                        high_EE_corr_diag{dir, trial}(ckpt) = prctile(diag(avg_EE), 100);
                        avg_EE_corr_subdiag{dir, trial}(ckpt) = mean(diag(avg_EE([2:corr_seg, 1], :)));
                        std_EE_corr_subdiag{dir, trial}(ckpt) = std(diag(avg_EE([2:corr_seg, 1], :)), 1);
                        low_EE_corr_subdiag{dir, trial}(ckpt) = prctile(diag(avg_EE([2:corr_seg, 1], :)), 0);
                        high_EE_corr_subdiag{dir, trial}(ckpt) = prctile(diag(avg_EE([2:corr_seg, 1], :)), 100);
                        avg_EE_uncorr{dir, trial}(ckpt) = (sum(sum(avg_EE)) - avg_EE_corr_diag{dir, trial}(ckpt)*corr_seg - ...
                                avg_EE_corr_subdiag{dir, trial}(ckpt)*corr_seg)/corr_seg/(corr_seg-2);
                        std_EE_uncorr{dir, trial}(ckpt) = sqrt((sum(sum(avg_EE.^2)) - sum(diag(avg_EE.^2)) - ...
                                sum(diag(avg_EE([2:corr_seg, 1], :).^2)))/corr_seg/(corr_seg-1) - ...
                                avg_EE_uncorr{dir, trial}(ckpt)^2);
                        low_EE_uncorr{dir, trial}(ckpt) = prctile(avg_EE(setdiff(1:corr_seg*corr_seg, [linspace(1, corr_seg*corr_seg, corr_seg), ...
                                linspace(2, corr_seg*(corr_seg-1), corr_seg-1), corr_seg*(corr_seg-1)+1])), 0);
                        high_EE_uncorr{dir, trial}(ckpt) = prctile(avg_EE(setdiff(1:corr_seg*corr_seg, [linspace(1, corr_seg*corr_seg, corr_seg), ...
                                linspace(2, corr_seg*(corr_seg-1), corr_seg-1), corr_seg*(corr_seg-1)+1])), 100);

                        avg_EE_all_reorder{dir, trial}(:, :, ckpt) = avg_EE_reorder;
                        avg_EE_corr_diag_reorder{dir, trial}(ckpt) = mean(diag(avg_EE_reorder));
                        std_EE_corr_diag_reorder{dir, trial}(ckpt) = std(diag(avg_EE_reorder), 1);
                        low_EE_corr_diag_reorder{dir, trial}(ckpt) = prctile(diag(avg_EE_reorder), 0);
                        high_EE_corr_diag_reorder{dir, trial}(ckpt) = prctile(diag(avg_EE_reorder), 100);
                        avg_EE_corr_subdiag_reorder{dir, trial}(ckpt) = mean(diag(avg_EE_reorder([2:corr_seg, 1], :)));
                        std_EE_corr_subdiag_reorder{dir, trial}(ckpt) = std(diag(avg_EE_reorder([2:corr_seg, 1], :)), 1);
                        low_EE_corr_subdiag_reorder{dir, trial}(ckpt) = prctile(diag(avg_EE_reorder([2:corr_seg, 1], :)), 0);
                        high_EE_corr_subdiag_reorder{dir, trial}(ckpt) = prctile(diag(avg_EE_reorder([2:corr_seg, 1], :)), 100);
                        avg_EE_uncorr_reorder{dir, trial}(ckpt) = (sum(sum(avg_EE_reorder)) - avg_EE_corr_diag_reorder{dir, trial}(ckpt)*corr_seg - ...
                                avg_EE_corr_subdiag_reorder{dir, trial}(ckpt)*corr_seg)/corr_seg/(corr_seg-2);
                        std_EE_uncorr_reorder{dir, trial}(ckpt) = sqrt((sum(sum(avg_EE_reorder.^2)) - sum(diag(avg_EE_reorder.^2)) - ...
                                sum(diag(avg_EE_reorder([2:corr_seg, 1], :).^2)))/corr_seg/(corr_seg-1) - ...
                                avg_EE_uncorr_reorder{dir, trial}(ckpt)^2);
                        low_EE_uncorr_reorder{dir, trial}(ckpt) = prctile(avg_EE_reorder(setdiff(1:corr_seg*corr_seg, [linspace(1, corr_seg*corr_seg, corr_seg), ...
                                linspace(2, corr_seg*(corr_seg-1), corr_seg-1), corr_seg*(corr_seg-1)+1])), 0);
                        high_EE_uncorr_reorder{dir, trial}(ckpt) = prctile(avg_EE_reorder(setdiff(1:corr_seg*corr_seg, [linspace(1, corr_seg*corr_seg, corr_seg), ...
                                linspace(2, corr_seg*(corr_seg-1), corr_seg-1), corr_seg*(corr_seg-1)+1])), 100);
                end
        end
end

save('trace_3/trace_3.mat', 'avg_EE_corr_diag', 'avg_EE_corr_subdiag', 'avg_EE_uncorr', ...
        'std_EE_corr_diag', 'std_EE_corr_subdiag', 'std_EE_uncorr', 'low_EE_corr_diag', ...
        'low_EE_corr_subdiag', 'low_EE_uncorr', 'high_EE_corr_diag', 'high_EE_corr_subdiag', ...
        'high_EE_uncorr', 'avg_EE_all', ...
        'avg_EE_corr_diag_reorder', 'avg_EE_corr_subdiag_reorder', 'avg_EE_uncorr_reorder', ...
        'std_EE_corr_diag_reorder', 'std_EE_corr_subdiag_reorder', 'std_EE_uncorr_reorder', 'low_EE_corr_diag_reorder', ...
        'low_EE_corr_subdiag_reorder', 'low_EE_uncorr_reorder', 'high_EE_corr_diag_reorder', 'high_EE_corr_subdiag_reorder', ...
        'high_EE_uncorr_reorder', 'avg_EE_all_reorder');
%}
load('trace_3/trace_3.mat');

avg_EE_connected = zeros(corr_seg, corr_seg, ckpt_num_all);
avg_EE_corr_diag_connected = zeros(ckpt_num_all, 1);
avg_EE_corr_subdiag_connected = zeros(ckpt_num_all, 1);
avg_EE_uncorr_connected = zeros(ckpt_num_all, 1);
std_EE_corr_diag_connected = zeros(ckpt_num_all, 1);
std_EE_corr_subdiag_connected = zeros(ckpt_num_all, 1);
std_EE_uncorr_connected = zeros(ckpt_num_all, 1);
low_EE_corr_diag_connected = zeros(ckpt_num_all, 1);
low_EE_corr_subdiag_connected = zeros(ckpt_num_all, 1);
low_EE_uncorr_connected = zeros(ckpt_num_all, 1);
high_EE_corr_diag_connected = zeros(ckpt_num_all, 1);
high_EE_corr_subdiag_connected = zeros(ckpt_num_all, 1);
high_EE_uncorr_connected = zeros(ckpt_num_all, 1);

avg_EE_connected_reorder = zeros(corr_seg, corr_seg, ckpt_num_all);
avg_EE_corr_diag_connected_reorder = zeros(ckpt_num_all, 1);
avg_EE_corr_subdiag_connected_reorder = zeros(ckpt_num_all, 1);
avg_EE_uncorr_connected_reorder = zeros(ckpt_num_all, 1);
std_EE_corr_diag_connected_reorder = zeros(ckpt_num_all, 1);
std_EE_corr_subdiag_connected_reorder = zeros(ckpt_num_all, 1);
std_EE_uncorr_connected_reorder = zeros(ckpt_num_all, 1);
low_EE_corr_diag_connected_reorder = zeros(ckpt_num_all, 1);
low_EE_corr_subdiag_connected_reorder = zeros(ckpt_num_all, 1);
low_EE_uncorr_connected_reorder = zeros(ckpt_num_all, 1);
high_EE_corr_diag_connected_reorder = zeros(ckpt_num_all, 1);
high_EE_corr_subdiag_connected_reorder = zeros(ckpt_num_all, 1);
high_EE_uncorr_connected_reorder = zeros(ckpt_num_all, 1);

start_frame = 1;

for dir = 1:length(dir_all)
        
        ckpt_num = round(duration_all(dir)/w_interval) + 1;
        
        avg_EE_dir = permute(reshape(cell2mat(avg_EE_all(dir, :)), [corr_seg, corr_seg, trials, ckpt_num]), [1, 2, 4, 3]);
        avg_EE_corr_diag_dir = cell2mat(avg_EE_corr_diag(dir, :));
        avg_EE_corr_subdiag_dir = cell2mat(avg_EE_corr_subdiag(dir, :));
        avg_EE_uncorr_dir = cell2mat(avg_EE_uncorr(dir, :));
        std_EE_corr_diag_dir = cell2mat(std_EE_corr_diag(dir, :));
        std_EE_corr_subdiag_dir = cell2mat(std_EE_corr_subdiag(dir, :));
        std_EE_uncorr_dir = cell2mat(std_EE_uncorr(dir, :));
        low_EE_corr_diag_dir = cell2mat(low_EE_corr_diag(dir, :));
        low_EE_corr_subdiag_dir = cell2mat(low_EE_corr_subdiag(dir, :));
        low_EE_uncorr_dir = cell2mat(low_EE_uncorr(dir, :));
        low_EE_uncorr_dir(isnan(low_EE_uncorr_dir)) = 0;
        high_EE_corr_diag_dir = cell2mat(high_EE_corr_diag(dir, :));
        high_EE_corr_subdiag_dir = cell2mat(high_EE_corr_subdiag(dir, :));
        high_EE_uncorr_dir = cell2mat(high_EE_uncorr(dir, :));
        high_EE_uncorr_dir(isnan(high_EE_uncorr_dir)) = 0;
        
        avg_EE_dir_reorder = permute(reshape(cell2mat(avg_EE_all_reorder(dir, :)), [corr_seg, corr_seg, trials, ckpt_num]), [1, 2, 4, 3]);
        avg_EE_corr_diag_dir_reorder = cell2mat(avg_EE_corr_diag_reorder(dir, :));
        avg_EE_corr_subdiag_dir_reorder = cell2mat(avg_EE_corr_subdiag_reorder(dir, :));
        avg_EE_uncorr_dir_reorder = cell2mat(avg_EE_uncorr_reorder(dir, :));
        std_EE_corr_diag_dir_reorder = cell2mat(std_EE_corr_diag_reorder(dir, :));
        std_EE_corr_subdiag_dir_reorder = cell2mat(std_EE_corr_subdiag_reorder(dir, :));
        std_EE_uncorr_dir_reorder = cell2mat(std_EE_uncorr_reorder(dir, :));
        low_EE_corr_diag_dir_reorder = cell2mat(low_EE_corr_diag_reorder(dir, :));
        low_EE_corr_subdiag_dir_reorder = cell2mat(low_EE_corr_subdiag_reorder(dir, :));
        low_EE_uncorr_dir_reorder = cell2mat(low_EE_uncorr_reorder(dir, :));
        low_EE_uncorr_dir_reorder(isnan(low_EE_uncorr_dir_reorder)) = 0;
        high_EE_corr_diag_dir_reorder = cell2mat(high_EE_corr_diag_reorder(dir, :));
        high_EE_corr_subdiag_dir_reorder = cell2mat(high_EE_corr_subdiag_reorder(dir, :));
        high_EE_uncorr_dir_reorder = cell2mat(high_EE_uncorr_reorder(dir, :));
        high_EE_uncorr_dir_reorder(isnan(high_EE_uncorr_dir_reorder)) = 0;
        
        if dir == 1
                avg_EE_connected(:, :, start_frame:(start_frame + size(avg_EE_all{dir}, 3) - 1)) = mean(avg_EE_dir, 4);
                avg_EE_corr_diag_connected(start_frame:(start_frame + size(avg_EE_all{dir}, 3) - 1)) = mean(avg_EE_corr_diag_dir, 2);
                avg_EE_corr_subdiag_connected(start_frame:(start_frame + size(avg_EE_all{dir}, 3) - 1)) = mean(avg_EE_corr_subdiag_dir, 2);
                avg_EE_uncorr_connected(start_frame:(start_frame + size(avg_EE_all{dir}, 3) - 1)) = mean(avg_EE_uncorr_dir, 2);
                std_EE_corr_diag_connected(start_frame:(start_frame + size(avg_EE_all{dir}, 3) - 1)) = sqrt(mean(std_EE_corr_diag_dir.^2, 2));
                std_EE_corr_subdiag_connected(start_frame:(start_frame + size(avg_EE_all{dir}, 3) - 1)) = sqrt(mean(std_EE_corr_subdiag_dir.^2, 2));
                std_EE_uncorr_connected(start_frame:(start_frame + size(avg_EE_all{dir}, 3) - 1)) = sqrt(mean(std_EE_uncorr_dir.^2, 2));
                low_EE_corr_diag_connected(start_frame:(start_frame + size(avg_EE_all{dir}, 3) - 1)) = mean(low_EE_corr_diag_dir, 2);
                low_EE_corr_subdiag_connected(start_frame:(start_frame + size(avg_EE_all{dir}, 3) - 1)) = mean(low_EE_corr_subdiag_dir, 2);
                low_EE_uncorr_connected(start_frame:(start_frame + size(avg_EE_all{dir}, 3) - 1)) = mean(low_EE_uncorr_dir, 2);
                high_EE_corr_diag_connected(start_frame:(start_frame + size(avg_EE_all{dir}, 3) - 1)) = mean(high_EE_corr_diag_dir, 2);
                high_EE_corr_subdiag_connected(start_frame:(start_frame + size(avg_EE_all{dir}, 3) - 1)) = mean(high_EE_corr_subdiag_dir, 2);
                high_EE_uncorr_connected(start_frame:(start_frame + size(avg_EE_all{dir}, 3) - 1)) = mean(high_EE_uncorr_dir, 2);
                
                avg_EE_connected_reorder(:, :, start_frame:(start_frame + size(avg_EE_all{dir}, 3) - 1)) = mean(avg_EE_dir_reorder, 4);
                avg_EE_corr_diag_connected_reorder(start_frame:(start_frame + size(avg_EE_all{dir}, 3) - 1)) = mean(avg_EE_corr_diag_dir_reorder, 2);
                avg_EE_corr_subdiag_connected_reorder(start_frame:(start_frame + size(avg_EE_all{dir}, 3) - 1)) = mean(avg_EE_corr_subdiag_dir_reorder, 2);
                avg_EE_uncorr_connected_reorder(start_frame:(start_frame + size(avg_EE_all{dir}, 3) - 1)) = mean(avg_EE_uncorr_dir_reorder, 2);
                std_EE_corr_diag_connected_reorder(start_frame:(start_frame + size(avg_EE_all{dir}, 3) - 1)) = sqrt(mean(std_EE_corr_diag_dir.^2, 2));
                std_EE_corr_subdiag_connected_reorder(start_frame:(start_frame + size(avg_EE_all{dir}, 3) - 1)) = sqrt(mean(std_EE_corr_subdiag_dir.^2, 2));
                std_EE_uncorr_connected_reorder(start_frame:(start_frame + size(avg_EE_all{dir}, 3) - 1)) = sqrt(mean(std_EE_uncorr_dir.^2, 2));
                low_EE_corr_diag_connected_reorder(start_frame:(start_frame + size(avg_EE_all{dir}, 3) - 1)) = mean(low_EE_corr_diag_dir_reorder, 2);
                low_EE_corr_subdiag_connected_reorder(start_frame:(start_frame + size(avg_EE_all{dir}, 3) - 1)) = mean(low_EE_corr_subdiag_dir_reorder, 2);
                low_EE_uncorr_connected_reorder(start_frame:(start_frame + size(avg_EE_all{dir}, 3) - 1)) = mean(low_EE_uncorr_dir_reorder, 2);
                high_EE_corr_diag_connected_reorder(start_frame:(start_frame + size(avg_EE_all{dir}, 3) - 1)) = mean(high_EE_corr_diag_dir_reorder, 2);
                high_EE_corr_subdiag_connected_reorder(start_frame:(start_frame + size(avg_EE_all{dir}, 3) - 1)) = mean(high_EE_corr_subdiag_dir_reorder, 2);
                high_EE_uncorr_connected_reorder(start_frame:(start_frame + size(avg_EE_all{dir}, 3) - 1)) = mean(high_EE_uncorr_dir_reorder, 2);
                
                start_frame = start_frame + size(avg_EE_all{dir}, 3);
        else
                avg_EE_connected(:, :, start_frame:(start_frame + size(avg_EE_all{dir}, 3) - 2)) = mean(avg_EE_dir(:, :, 2:end, :), 4);
                avg_EE_corr_diag_connected(start_frame:(start_frame + size(avg_EE_all{dir}, 3) - 2)) = mean(avg_EE_corr_diag_dir(2:end, :), 2);
                avg_EE_corr_subdiag_connected(start_frame:(start_frame + size(avg_EE_all{dir}, 3) - 2)) = mean(avg_EE_corr_subdiag_dir(2:end, :), 2);
                avg_EE_uncorr_connected(start_frame:(start_frame + size(avg_EE_all{dir}, 3) - 2)) = mean(avg_EE_uncorr_dir(2:end, :), 2);
                std_EE_corr_diag_connected(start_frame:(start_frame + size(avg_EE_all{dir}, 3) - 2)) = sqrt(mean(std_EE_corr_diag_dir(2:end, :).^2, 2));
                std_EE_corr_subdiag_connected(start_frame:(start_frame + size(avg_EE_all{dir}, 3) - 2)) = sqrt(mean(std_EE_corr_subdiag_dir(2:end, :).^2, 2));
                std_EE_uncorr_connected(start_frame:(start_frame + size(avg_EE_all{dir}, 3) - 2)) = sqrt(mean(std_EE_uncorr_dir(2:end, :).^2, 2));
                low_EE_corr_diag_connected(start_frame:(start_frame + size(avg_EE_all{dir}, 3) - 2)) = mean(low_EE_corr_diag_dir(2:end, :), 2);
                low_EE_corr_subdiag_connected(start_frame:(start_frame + size(avg_EE_all{dir}, 3) - 2)) = mean(low_EE_corr_subdiag_dir(2:end, :), 2);
                low_EE_uncorr_connected(start_frame:(start_frame + size(avg_EE_all{dir}, 3) - 2)) = mean(low_EE_uncorr_dir(2:end, :), 2);
                high_EE_corr_diag_connected(start_frame:(start_frame + size(avg_EE_all{dir}, 3) - 2)) = mean(high_EE_corr_diag_dir(2:end, :), 2);
                high_EE_corr_subdiag_connected(start_frame:(start_frame + size(avg_EE_all{dir}, 3) - 2)) = mean(high_EE_corr_subdiag_dir(2:end, :), 2);
                high_EE_uncorr_connected(start_frame:(start_frame + size(avg_EE_all{dir}, 3) - 2)) = mean(high_EE_uncorr_dir(2:end, :), 2);
                
                avg_EE_connected_reorder(:, :, start_frame:(start_frame + size(avg_EE_all{dir}, 3) - 2)) = mean(avg_EE_dir_reorder(:, :, 2:end, :), 4);
                avg_EE_corr_diag_connected_reorder(start_frame:(start_frame + size(avg_EE_all{dir}, 3) - 2)) = mean(avg_EE_corr_diag_dir_reorder(2:end, :), 2);
                avg_EE_corr_subdiag_connected_reorder(start_frame:(start_frame + size(avg_EE_all{dir}, 3) - 2)) = mean(avg_EE_corr_subdiag_dir_reorder(2:end, :), 2);
                avg_EE_uncorr_connected_reorder(start_frame:(start_frame + size(avg_EE_all{dir}, 3) - 2)) = mean(avg_EE_uncorr_dir_reorder(2:end, :), 2);
                std_EE_corr_diag_connected_reorder(start_frame:(start_frame + size(avg_EE_all{dir}, 3) - 2)) = sqrt(mean(std_EE_corr_diag_dir_reorder(2:end, :).^2, 2));
                std_EE_corr_subdiag_connected_reorder(start_frame:(start_frame + size(avg_EE_all{dir}, 3) - 2)) = sqrt(mean(std_EE_corr_subdiag_dir_reorder(2:end, :).^2, 2));
                std_EE_uncorr_connected_reorder(start_frame:(start_frame + size(avg_EE_all{dir}, 3) - 2)) = sqrt(mean(std_EE_uncorr_dir_reorder(2:end, :).^2, 2));
                low_EE_corr_diag_connected_reorder(start_frame:(start_frame + size(avg_EE_all{dir}, 3) - 2)) = mean(low_EE_corr_diag_dir_reorder(2:end, :), 2);
                low_EE_corr_subdiag_connected_reorder(start_frame:(start_frame + size(avg_EE_all{dir}, 3) - 2)) = mean(low_EE_corr_subdiag_dir_reorder(2:end, :), 2);
                low_EE_uncorr_connected_reorder(start_frame:(start_frame + size(avg_EE_all{dir}, 3) - 2)) = mean(low_EE_uncorr_dir_reorder(2:end, :), 2);
                high_EE_corr_diag_connected_reorder(start_frame:(start_frame + size(avg_EE_all{dir}, 3) - 2)) = mean(high_EE_corr_diag_dir_reorder(2:end, :), 2);
                high_EE_corr_subdiag_connected_reorder(start_frame:(start_frame + size(avg_EE_all{dir}, 3) - 2)) = mean(high_EE_corr_subdiag_dir_reorder(2:end, :), 2);
                high_EE_uncorr_connected_reorder(start_frame:(start_frame + size(avg_EE_all{dir}, 3) - 2)) = mean(high_EE_uncorr_dir_reorder(2:end, :), 2);
                
                start_frame = start_frame + size(avg_EE_all{dir}, 3) - 1;
        end
end



figure
plot((0:(ckpt_num_all-1))*w_interval, avg_EE_corr_diag_connected, 'c', 'LineWidth', 2)
hold on
plot((0:(ckpt_num_all-1))*w_interval, avg_EE_corr_subdiag_connected, 'm', 'LineWidth', 2)
plot((0:(ckpt_num_all-1))*w_interval, avg_EE_uncorr_connected, 'color', [1 1 1]*0.3, 'LineWidth', 2)
plot((0:(ckpt_num_all-1))*w_interval, avg_EE_corr_diag_connected_reorder, 'b', 'LineWidth', 2)
plot((0:(ckpt_num_all-1))*w_interval, avg_EE_corr_subdiag_connected_reorder, 'r', 'LineWidth', 2)
plot((0:(ckpt_num_all-1))*w_interval, avg_EE_uncorr_connected_reorder, 'k', 'LineWidth', 2)
fill([0:(ckpt_num_all-1), flip(0:(ckpt_num_all-1))]*w_interval, [low_EE_corr_diag_connected; ...
        flip(high_EE_corr_diag_connected)], [1 1 1]*0.8, 'LineStyle', 'none')
fill([0:(ckpt_num_all-1), flip(0:(ckpt_num_all-1))]*w_interval, [low_EE_corr_subdiag_connected; ...
        flip(high_EE_corr_subdiag_connected)], [1 1 1]*0.8, 'LineStyle', 'none')
fill([0:(ckpt_num_all-1), flip(0:(ckpt_num_all-1))]*w_interval, [low_EE_corr_diag_connected_reorder; ...
        flip(high_EE_corr_diag_connected_reorder)], [0.5 0.5 1], 'LineStyle', 'none')
fill([0:(ckpt_num_all-1), flip(0:(ckpt_num_all-1))]*w_interval, [low_EE_corr_subdiag_connected_reorder; ...
        flip(high_EE_corr_subdiag_connected_reorder)], [1 0.5 0.5], 'LineStyle', 'none')
fill([0:(ckpt_num_all-1), flip(0:(ckpt_num_all-1))]*w_interval, [low_EE_uncorr_connected; ...
        flip(high_EE_uncorr_connected)], [1 1 1]*0.8, 'LineStyle', 'none')
fill([0:(ckpt_num_all-1), flip(0:(ckpt_num_all-1))]*w_interval, [low_EE_uncorr_connected_reorder; ...
        flip(high_EE_uncorr_connected_reorder)], [1 1 1]*0.8, 'LineStyle', 'none')
plot((0:(ckpt_num_all-1))*w_interval, avg_EE_corr_diag_connected, 'c', 'LineWidth', 2)
plot((0:(ckpt_num_all-1))*w_interval, avg_EE_corr_subdiag_connected, 'm', 'LineWidth', 2)
plot((0:(ckpt_num_all-1))*w_interval, avg_EE_corr_diag_connected_reorder, 'b', 'LineWidth', 2)
plot((0:(ckpt_num_all-1))*w_interval, avg_EE_corr_subdiag_connected_reorder, 'r', 'LineWidth', 2)
plot((0:(ckpt_num_all-1))*w_interval, avg_EE_uncorr_connected, 'color', [1 1 1]*0.3, 'LineWidth', 2)
plot((0:(ckpt_num_all-1))*w_interval, avg_EE_uncorr_connected_reorder, 'k', 'LineWidth', 2)
set(gcf, 'color', [1 1 1])   
ax = gca;
a = get(ax,'Label');
set(gca,'Label', a,'fontsize',15)     
xlabel('Time (s)')
ylabel('Ratio of high effective weight');
legend('Old diagonal', 'Old subdiagonal', 'Old uncorrelated', 'New diagonal', 'New subdiagonal', 'New uncorrelated', 'Location', 'northeast')
savefig(gcf, 'trace_3/avg_wth_ratio_across_time.fig')

figure
plot((0:(ckpt_num_all-1))*w_interval, avg_EE_corr_diag_connected, 'b', 'LineWidth', 2)
hold on
fill([0:(ckpt_num_all-1), flip(0:(ckpt_num_all-1))]*w_interval, [low_EE_corr_diag_connected; ...
        flip(high_EE_corr_diag_connected)]', [0.7 0.7 1], 'LineStyle', 'none')
plot((0:(ckpt_num_all-1))*w_interval, avg_EE_corr_subdiag_connected, 'r', 'LineWidth', 2)
fill([0:(ckpt_num_all-1), flip(0:(ckpt_num_all-1))]*w_interval, [low_EE_corr_subdiag_connected; ...
        flip(high_EE_corr_subdiag_connected)]', [1 0.7 0.7], 'LineStyle', 'none')
plot((0:(ckpt_num_all-1))*w_interval, avg_EE_uncorr_connected, 'k', 'LineWidth', 2)
fill([0:(ckpt_num_all-1), flip(0:(ckpt_num_all-1))]*w_interval, [low_EE_uncorr_connected; ...
        flip(high_EE_uncorr_connected)]', [0.7 0.7 0.7], 'LineStyle', 'none')
plot((0:(ckpt_num_all-1))*w_interval, avg_EE_corr_diag_connected, 'b', 'LineWidth', 2)
hold on
plot((0:(ckpt_num_all-1))*w_interval, avg_EE_corr_subdiag_connected, 'r', 'LineWidth', 2)
plot((0:(ckpt_num_all-1))*w_interval, avg_EE_uncorr_connected, 'k', 'LineWidth', 2)
set(gcf, 'color', [1 1 1])   
ax = gca;
a = get(ax,'Label');
set(gca,'Label', a,'fontsize',15)     
xlabel('Time (s)')
ylabel('Ratio of high effective weight');
legend('Diagonal', 'Whole range', 'Subdiagonal', 'Whole range', 'Uncorrelated', 'Whole range', 'Location', 'northeast')
savefig(gcf, 'trace_3/avg_wth_ratio_across_time_2.fig')

start_ckpt = 1;
ckpt_reorder = zeros(ckpt_num_all, 1);
for dir = 1:length(dir_all)
        if dir == 1
                ckpt_reorder(1:(round(duration_all(dir)/w_interval)+1)) = reorder(dir);
        else
                ckpt_reorder((round(sum(duration_all(1:(dir-1)))/w_interval)+2): round(sum(duration_all(1:dir))/w_interval)+1) = reorder(dir);
        end
end

for ckpt = 1:ckpt_num_all
        if mod(ckpt, 10) == 1
                figure
                if ckpt_reorder(ckpt) == 1
                        imagesc(avg_EE_connected_reorder(:, :, ckpt));
                else
                        imagesc(avg_EE_connected(:, :, ckpt));
                end
                set(gcf, 'color', [1 1 1])
                xticks([4, 8, 12])
                yticks([4, 8, 12]);
                xticklabels({'400', '800', '1200'})
                yticklabels({'400', '800', '1200'})
                ax = gca;
                a = get(ax,'Label');
                set(gca,'Label', a,'fontsize',15)
                xlabel('#Neuron')
                ylabel('#Neuron')
                h = colorbar;
                set(get(h, 'label'), 'string', 'Ratio of high effective weight');
                %h.Limits = [0 7e-3];
                title(['t = ', num2str(w_interval*(ckpt-1)), 's']);
                savefig(gcf, ['trace_3/avg_wth_ratio_', num2str(w_interval*(ckpt-1)), 's']);
        end
end

save('trace_3/trace_3_plot.mat', 'ckpt_num_all', 'w_interval', 'avg_EE_corr_diag_connected', 'avg_EE_corr_subdiag_connected', ...
        'avg_EE_uncorr_connected', 'low_EE_corr_diag_connected', 'high_EE_corr_diag_connected', 'low_EE_corr_subdiag_connected', ...
        'high_EE_corr_subdiag_connected', 'low_EE_uncorr_connected', 'high_EE_uncorr_connected', ...
        'avg_EE_corr_diag_connected', 'avg_EE_corr_subdiag_connected', 'avg_EE_uncorr_connected', 'avg_EE_connected', ...
        'avg_EE_corr_diag_connected_reorder', 'avg_EE_corr_subdiag_connected_reorder', 'avg_EE_uncorr_connected_reorder', ...
        'low_EE_corr_diag_connected_reorder', 'high_EE_corr_diag_connected_reorder', 'low_EE_corr_subdiag_connected_reorder', ...
        'high_EE_corr_subdiag_connected_reorder', 'low_EE_uncorr_connected_reorder', 'high_EE_uncorr_connected_reorder', ...
        'avg_EE_corr_diag_connected_reorder', 'avg_EE_corr_subdiag_connected_reorder', ...
        'avg_EE_uncorr_connected_reorder', 'avg_EE_connected_reorder')
                
                
%}
        