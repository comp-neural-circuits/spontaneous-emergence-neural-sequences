%samples = [584, 406, 1129, 387, 153, 237, 366, 715, 738, 969];
%sample_groups = [22, 21, 19, 19, 18, 18, 18, 18, 18, 18];

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

figure
set(gcf, 'color', [1 1 1])

norm_entropy_data_all = - ones(length(samples), cutoff)*NaN;
norm_entropy_shuffle_all = - ones(length(samples)*shuffle_num, cutoff)*NaN;

for i = 1:length(samples)
        load([num2str(sample_groups(i)), '/', num2str(samples(i)), '/sequence/entropy.mat']);
        norm_entropy_data_all(i, 1:min(cutoff, length(norm_entropy_data))) = ...
                norm_entropy_data(1:min(cutoff, length(norm_entropy_data)));
        norm_entropy_shuffle_all((i-1)*shuffle_num+(1:shuffle_num), ...
                1:min(cutoff, length(norm_entropy_data))) = ...
                norm_entropy_shuffle(:, 1:min(cutoff, length(norm_entropy_data)));
        %{
        plot(1:min(cutoff, length(norm_entropy_data)), norm_entropy_data(...
                1:min(cutoff, length(norm_entropy_data))), 'k', 'LineWidth', 1);
        hold on
        plot(1:min(cutoff, length(norm_entropy_data)), norm_entropy_shuffle(...
                :, 1:min(cutoff, length(norm_entropy_data))), 'r', 'LineWidth', 1)
                %}
end

norm_entropy_data_avg = zeros(1, cutoff);
norm_entropy_shuffle_avg = zeros(1, cutoff);
norm_entropy_data_std = zeros(1, cutoff);
norm_entropy_shuffle_std = zeros(1, cutoff);

for i = 1:cutoff
        norm_entropy_data_avg(i) = nanmean(norm_entropy_data_all(:, i));
        norm_entropy_shuffle_avg(i) = nanmean(norm_entropy_shuffle_all(:, i));
        norm_entropy_data_std(i) = nanstd(norm_entropy_data_all(:, i));
        norm_entropy_shuffle_std(i) = nanstd(norm_entropy_shuffle_all(:, i));
end

save('entropy_rank_correlation/entropy.mat', 'cutoff', 'norm_entropy_data_avg', ...
        'norm_entropy_shuffle_avg', 'norm_entropy_data_std', 'norm_entropy_shuffle_std')

%plot(norm_entropy_data_avg, 'k', 'LineWidth', 3)
%plot(norm_entropy_shuffle_avg, 'r', 'LineWidth', 3)
fill([(1:cutoff), flip(1:cutoff)]', [norm_entropy_data_avg - norm_entropy_data_std, ...
        flip(norm_entropy_data_avg + norm_entropy_data_std)]', [1 1 1]*0.7, 'LineStyle', 'none')
hold on
fill([(1:cutoff), flip(1:cutoff)]', [norm_entropy_shuffle_avg - norm_entropy_shuffle_std, ...
        flip(norm_entropy_shuffle_avg + norm_entropy_shuffle_std)]', [1 0.7 0.7], 'LineStyle', 'none')
plot(norm_entropy_data_avg, 'k', 'LineWidth', 3)
plot(norm_entropy_shuffle_avg, 'r', 'LineWidth', 3)
xlabel('Spike rank')
ylabel('Norm. entropy')
%legend('Simulations', 'Shuffled', 'Location', 'southeast')
legend('Std', 'Std', 'Model', 'Shuffled', 'Location', 'southeast')
ax = gca;
a = get(ax,'Label');
set(gca,'Label', a,'fontsize',15)
ylim([0 1])
set(gca, 'Box', 'off')
        
savefig(gcf, 'entropy_rank_correlation/entropy.fig')
        