groups = 15;

followers_corr = cell(groups, 1);
followers_uncorr = cell(groups, 1);
avg_prob_corr = cell(groups, 1);
avg_prob_uncorr = cell(groups, 1);
all_prob_corr = cell(groups, 1);
all_prob_uncorr = cell(groups, 1);
avg_delay_corr = cell(groups, 1);
avg_delay_uncorr = cell(groups, 1);
all_delay_corr = cell(groups, 1);
all_delay_uncorr = cell(groups, 1);

corr_range = 100;

for i = 1:groups
        cd(num2str(i))
        folders = dir();
        folders = folders(3:end);
        followers_corr{i} = zeros(length(folders), 1);
        followers_uncorr{i} = zeros(length(folders), 1);
        for j = 1:length(folders)
                cd(folders(j).name)
                load('sequence/delay_and_jitter.mat')
                source = str2double(folders(j).name);
                corr_seg = ceil(source/corr_range);
                followers_corr{i}(j) = sum(index_e >= (1+(corr_seg - 1)*corr_range) & index_e <= corr_seg*corr_range) - 1;
                followers_uncorr{i}(j) = length(index_e) - 1 - followers_corr{i}(j);
                idx_corr = index_e >= (1+(corr_seg - 1)*corr_range) & index_e <= corr_seg*corr_range & index_e ~= source;
                idx_uncorr = index_e <= (corr_seg - 1)*corr_range | index_e > corr_seg*corr_range;
                avg_prob_corr{i}(j) = mean(prob_follower_e(idx_corr));
                avg_prob_uncorr{i}(j) = mean(prob_follower_e(idx_uncorr));
                all_prob_corr{i} = [all_prob_corr{i}; prob_follower_e(idx_corr)];
                all_prob_uncorr{i} = [all_prob_uncorr{i}; prob_follower_e(idx_uncorr)];
                avg_delay_corr{i}(j) = mean(med_delay_e(idx_corr));
                avg_delay_uncorr{i}(j) = mean(med_delay_e(idx_uncorr));
                all_delay_corr{i} = [all_delay_corr{i}; med_delay_e(idx_corr)];
                all_delay_uncorr{i} = [all_delay_uncorr{i}; med_delay_e(idx_uncorr)];
                cd('..');
        end
        cd('..')
end

save('followers_corr_uncorr.mat', 'followers_corr', 'followers_uncorr', 'avg_prob_corr', 'avg_prob_uncorr', ...
        'all_prob_corr', 'all_prob_uncorr', 'avg_delay_corr', 'avg_delay_uncorr', 'all_delay_corr', 'all_delay_uncorr');


