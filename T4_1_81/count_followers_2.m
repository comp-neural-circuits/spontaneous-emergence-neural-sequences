groups = 19;

followers_corr_diag = cell(groups, 1);
followers_corr_subdiag = cell(groups, 1);
followers_uncorr = cell(groups, 1);
avg_prob_corr_diag = cell(groups, 1);
avg_prob_corr_subdiag = cell(groups, 1);
avg_prob_uncorr = cell(groups, 1);
all_prob_corr_diag = cell(groups, 1);
all_prob_corr_subdiag = cell(groups, 1);
all_prob_uncorr = cell(groups, 1);
avg_delay_corr_diag = cell(groups, 1);
avg_delay_corr_subdiag = cell(groups, 1);
avg_delay_uncorr = cell(groups, 1);
all_delay_corr_diag = cell(groups, 1);
all_delay_corr_subdiag = cell(groups, 1);
all_delay_uncorr = cell(groups, 1);

corr_range = 100;

for i = 1:groups
        cd(num2str(i))
        folders = dir();
        folders = folders(3:end);
        followers_corr_diag{i} = zeros(length(folders), 1);
        followers_corr_subdiag{i} = zeros(length(folders), 1);
        followers_uncorr{i} = zeros(length(folders), 1);
        for j = 1:length(folders)
                cd(folders(j).name)
                load('sequence/delay_and_jitter.mat')
                source = str2double(folders(j).name);
                corr_seg = ceil(source/corr_range);
                followers_corr_diag{i}(j) = sum(index_e >= (1+(corr_seg - 1)*corr_range) & index_e <= corr_seg*corr_range) - 1;
                if corr_seg == 12
                        followers_corr_subdiag{i}(j) = sum(index_e <= corr_range);
                else 
                        followers_corr_subdiag{i}(j) = sum(index_e >= (1+corr_seg*corr_range) & index_e <= (corr_seg+1)*corr_range);
                end
                followers_uncorr{i}(j) = length(index_e) - 1 - followers_corr_diag{i}(j) - followers_corr_subdiag{i}(j);
                idx_corr_diag = index_e >= (1+(corr_seg - 1)*corr_range) & index_e <= corr_seg*corr_range & index_e ~= source;
                if corr_seg == 12
                        idx_corr_subdiag = index_e <= corr_range;
                        idx_uncorr = index_e <= (corr_seg - 1)*corr_range & index_e > corr_range;
                else
                        idx_corr_subdiag = index_e >= (1+corr_seg*corr_range) & index_e <= (corr_seg+1)*corr_range;
                        idx_uncorr = index_e <= (corr_seg - 1)*corr_range | index_e > (corr_seg+1)*corr_range;
                end
                avg_prob_corr_diag{i}(j) = mean(prob_follower_e(idx_corr_diag));
                avg_prob_corr_subdiag{i}(j) = mean(prob_follower_e(idx_corr_subdiag));
                avg_prob_uncorr{i}(j) = mean(prob_follower_e(idx_uncorr));
                all_prob_corr_diag{i} = [all_prob_corr_diag{i}; prob_follower_e(idx_corr_diag)];
                all_prob_corr_subdiag{i} = [all_prob_corr_subdiag{i}; prob_follower_e(idx_corr_subdiag)];
                all_prob_uncorr{i} = [all_prob_uncorr{i}; prob_follower_e(idx_uncorr)];
                avg_delay_corr_diag{i}(j) = mean(med_delay_e(idx_corr_diag));
                avg_delay_corr_subdiag{i}(j) = mean(med_delay_e(idx_corr_subdiag));
                avg_delay_uncorr{i}(j) = mean(med_delay_e(idx_uncorr));
                all_delay_corr_diag{i} = [all_delay_corr_diag{i}; med_delay_e(idx_corr_diag)];
                all_delay_corr_subdiag{i} = [all_delay_corr_subdiag{i}; med_delay_e(idx_corr_subdiag)];
                all_delay_uncorr{i} = [all_delay_uncorr{i}; med_delay_e(idx_uncorr)];
                cd('..');
        end
        cd('..')
end

save('followers_corr_uncorr2.mat', 'followers_corr_diag', 'followers_corr_subdiag', 'followers_uncorr', ...
        'avg_prob_corr_diag', 'avg_prob_corr_subdiag', 'avg_prob_uncorr', ...
        'all_prob_corr_diag', 'all_prob_corr_subdiag', 'all_prob_uncorr', ...
        'avg_delay_corr_diag', 'avg_delay_corr_subdiag', 'avg_delay_uncorr', ...
        'all_delay_corr_diag', 'all_delay_corr_subdiag', 'all_delay_uncorr');


