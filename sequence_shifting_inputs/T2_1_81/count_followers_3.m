%count followers in every block

groups = 19;
corr_range = 100;
Ne = 1200;
seg_num = round(Ne/corr_range);

num_followers = cell(groups, seg_num);
avg_prob = cell(groups, seg_num);
all_prob = cell(groups, seg_num);
avg_delay = cell(groups, seg_num);
all_delay = cell(groups, seg_num);

for i = 1:groups
        cd(num2str(i))
        folders = dir();
        folders = folders(3:end);
        for seg = 1:seg_num
                num_followers{i, seg} = zeros(length(folders), 1);
        end
        for j = 1:length(folders)
                cd(folders(j).name)
                load('sequence/delay_and_jitter.mat')
                source = str2double(folders(j).name);
                seg_source = ceil(source/corr_range);
                for seg = 1:seg_num
                        idx_start = 1+ mod((seg_source + seg - 2)*corr_range, Ne);
                        idx_end = idx_start + corr_range - 1;
                        idx_seg = index_e >= idx_start & index_e <= idx_end & index_e ~= source;
                        num_followers{i, seg}(j) = sum(idx_seg);
                        avg_prob{i, seg}(j) = mean(prob_follower_e(idx_seg));
                        all_prob{i, seg} = [all_prob{i, seg}; prob_follower_e(idx_seg)];
                        avg_delay{i, seg}(j) = mean(med_delay_e(idx_seg));
                        all_delay{i, seg} = [all_delay{i, seg}; med_delay_e(idx_seg)];
                end
                cd('..');
        end
        cd('..')
end

save('followers_all_blocks.mat', 'num_followers', 'avg_prob', 'all_prob', 'avg_delay', 'all_delay')


