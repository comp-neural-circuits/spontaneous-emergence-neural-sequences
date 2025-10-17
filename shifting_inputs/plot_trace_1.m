%6(5000), 7(3000), 10(7000)

dir_all = {'1', '2'};
duration_all = [80, 80]*1e2;
t_start = [0, 100];

p_ee = 0.17/3;

trials = 10;

w_interval = 100;

V0 = -0.063;

th_0 = 2.5;  %threshold for high weight-thresh ratio

Ne = 1200;

corr_range = 100;

corr_seg = round(Ne/corr_range);

ckpt_num_all = round(sum(duration_all)/w_interval)+1;

if ~isfolder('trace_1')
        mkdir('trace_1')
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

for dir = 1:length(dir_all)
        parfor trial = 1:trials
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
                        for i = 1:corr_seg
                                for j = 1:corr_seg
                                        avg_EE(i, j) = mean(mean(effw((1:corr_range)+(i-1)*corr_range, ...
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
                end
        end
end

save('trace_1/trace_1.mat', 'avg_EE_corr_diag', 'avg_EE_corr_subdiag', 'avg_EE_uncorr', ...
        'std_EE_corr_diag', 'std_EE_corr_subdiag', 'std_EE_uncorr', 'low_EE_corr_diag', ...
        'low_EE_corr_subdiag', 'low_EE_uncorr', 'high_EE_corr_diag', 'high_EE_corr_subdiag', ...
        'high_EE_uncorr', 'avg_EE_all');
%}
load('trace_1/trace_1.mat');

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
        high_EE_corr_diag_dir = cell2mat(high_EE_corr_diag(dir, :));
        high_EE_corr_subdiag_dir = cell2mat(high_EE_corr_subdiag(dir, :));
        high_EE_uncorr_dir = cell2mat(high_EE_uncorr(dir, :));
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
                start_frame = start_frame + size(avg_EE_all{dir}, 3) - 1;
        end
end

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
xlim([0 1.6e4])
legend('Diagonal', 'Whole range', 'Subdiagonal', 'Whole range', 'Uncorrelated', 'Whole range', 'Location', 'northeast')
savefig(gcf, 'trace_1/avg_wth_ratio_across_time.fig')

for ckpt = 1:ckpt_num_all
        if mod(ckpt, 20) == 1
                figure
                imagesc(avg_EE_connected(:, :, ckpt));
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
                h.Limits = [0 0.06];
                title(['t = ', num2str(w_interval*(ckpt-1)), 's']);
                savefig(gcf, ['trace_1/avg_wth_ratio_', num2str(w_interval*(ckpt-1)), 's']);
        end
end

save('trace_1/trace_1_plot.mat', 'ckpt_num_all', 'w_interval', 'avg_EE_corr_diag_connected', 'avg_EE_corr_subdiag_connected', ...
        'avg_EE_uncorr_connected', 'low_EE_corr_diag_connected', 'high_EE_corr_diag_connected', 'low_EE_corr_subdiag_connected', ...
        'high_EE_corr_subdiag_connected', 'low_EE_uncorr_connected', 'high_EE_uncorr_connected', ...
        'avg_EE_corr_diag_connected', 'avg_EE_corr_subdiag_connected', 'avg_EE_uncorr_connected', 'avg_EE_connected')
                
                
%}
        