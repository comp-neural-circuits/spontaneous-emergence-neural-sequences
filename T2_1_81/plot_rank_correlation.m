load('followers.mat')

if ~isfolder('entropy_rank_correlation')
        mkdir('entropy_rank_correlation')
end

folders = dir();

samples = [];
sample_groups = [];
for id = 3:(2+length(followers_e))
        cd(folders(id).name)
        sub_folders = dir();
        sub_folders = sub_folders(3:end);
        for sub_id = 1:length(sub_folders)
                if followers_e{str2double(folders(id).name)}(sub_id) + followers_i{str2double(folders(id).name)}(sub_id) > 0
                        samples = [samples, str2double(sub_folders(sub_id).name)];
                        sample_groups = [sample_groups, str2double(folders(id).name)];
                end
        end
        cd('..')
end

cutoff = 50;
shuffle_num = 10;

rank_corr_data_all = zeros(length(samples), cutoff)*NaN;
rank_corr_shuffle_all = zeros(length(samples)*shuffle_num, cutoff)*NaN;

for i = 1:length(samples)
        load([num2str(sample_groups(i)), '/', num2str(samples(i)), '/sequence/rank_corr.mat']);
        r_data_sort = r_data(sorted_neuron_idx);
        r_shuffle_sort = r_shuffle(:, sorted_neuron_idx);
        rank_corr_data_all(i, 1:min(cutoff, length(r_data))) = r_data_sort(1:min(cutoff, length(r_data)));
        rank_corr_shuffle_all((i-1)*shuffle_num+(1:shuffle_num), ...
                1:min(cutoff, length(r_data))) = r_shuffle_sort(:, 1:min(cutoff, length(r_data)));
        %plot(1:min(cutoff, length(r_data)), r_data_sort(1:min(cutoff, length(r_data))), 'color', [1 1 1]*0.5, 'LineWidth', 1);
        %hold on
        %plot(1:min(cutoff, length(r_data)), r_shuffle_sort(:, 1:min(cutoff, length(r_data))), 'color', [1 0.5 0.5],'LineWidth', 1)
end

rank_corr_data_avg = zeros(1, cutoff);
rank_corr_shuffle_avg = zeros(1, cutoff);
rank_corr_data_std = zeros(1, cutoff);
rank_corr_shuffle_std = zeros(1, cutoff);

for i = 1:cutoff
        rank_corr_data_avg(i) = nanmean(rank_corr_data_all(:, i));
        rank_corr_shuffle_avg(i) = nanmean(rank_corr_shuffle_all(:, i));
        rank_corr_data_std(i) = nanstd(rank_corr_data_all(:, i));
        rank_corr_shuffle_std(i) = nanstd(rank_corr_shuffle_all(:, i));
end

save('entropy_rank_correlation/rank_corr.mat', 'cutoff', 'rank_corr_data_avg', ...
        'rank_corr_shuffle_avg', 'rank_corr_data_std', 'rank_corr_shuffle_std')

figure
fill([(1:cutoff), flip(1:cutoff)]', [rank_corr_data_avg - rank_corr_data_std, ...
        flip(rank_corr_data_avg + rank_corr_data_std)]', [1 1 1]*0.7, 'LineStyle', 'none')
hold on
fill([(1:cutoff), flip(1:cutoff)]', [rank_corr_shuffle_avg - rank_corr_shuffle_std, ...
        flip(rank_corr_shuffle_avg + rank_corr_shuffle_std)]', [1 0.7 0.7], 'LineStyle', 'none')
plot(rank_corr_data_avg, 'k', 'LineWidth', 3)
plot(rank_corr_shuffle_avg, 'r', 'LineWidth', 3)
xlim([1 50])
xlabel('Spike rank')
ylabel('Rank Corr.')
legend('Std', 'Std', 'Simulation', 'Shuffled')
set(gcf, 'color', [1 1 1])
ax = gca;
a = get(ax,'Label');
set(gca,'Label', a,'fontsize',15)
set(gca, 'Box', 'off')
savefig(gcf, 'entropy_rank_correlation/rank_corr.fig');



