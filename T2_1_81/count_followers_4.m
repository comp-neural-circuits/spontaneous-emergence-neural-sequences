%count second next block

groups = 19;

followers_corr_diag = cell(groups, 1);
followers_corr_subdiag = cell(groups, 1);
followers_corr_subdiag2 = cell(groups, 1);
followers_uncorr = cell(groups, 1);
avg_prob_corr_diag = cell(groups, 1);
avg_prob_corr_subdiag = cell(groups, 1);
avg_prob_corr_subdiag2 = cell(groups, 1);
avg_prob_uncorr = cell(groups, 1);
all_prob_corr_diag = cell(groups, 1);
all_prob_corr_subdiag = cell(groups, 1);
all_prob_corr_subdiag2 = cell(groups, 1);
all_prob_uncorr = cell(groups, 1);
avg_delay_corr_diag = cell(groups, 1);
avg_delay_corr_subdiag = cell(groups, 1);
avg_delay_corr_subdiag2 = cell(groups, 1);
avg_delay_uncorr = cell(groups, 1);
all_delay_corr_diag = cell(groups, 1);
all_delay_corr_subdiag = cell(groups, 1);
all_delay_corr_subdiag2 = cell(groups, 1);
all_delay_uncorr = cell(groups, 1);

corr_range = 100;

Ne = 1200;

EE_dir = '../shifting_inputs/2/1';
EE_seg = h5read_sparsematrix([EE_dir, '/all_w.h5'], [num2str(81), '/EE']);
ET = h5read([EE_dir, '/all_w.h5'], ['/', num2str(81), '/ET']);
EE = retrieve_sp(EE_seg);
V0 = -0.063;
effw = EE./((ET - V0)*ones(1,size(EE, 2)))/1e3;

for i = 1:groups
        cd(num2str(i))
        folders = dir();
        folders = folders(3:end);
        followers_corr_diag{i} = zeros(length(folders), 1);
        followers_corr_subdiag{i} = zeros(length(folders), 1);
        followers_corr_subdiag2{i} = zeros(length(folders), 1);
        followers_uncorr{i} = zeros(length(folders), 1);
        for j = 1:length(folders)
                cd(folders(j).name)
                load('sequence/delay_and_jitter.mat')
                source = str2double(folders(j).name);
                corr_seg = ceil(source/corr_range);
                %index_fol_corr_diag = setdiff(index(index_e >= (1+(corr_seg - 1)*corr_range) & index_e <= corr_seg*corr_range), source);
                idx_corr_diag = find(index_e >= (1+(corr_seg - 1)*corr_range) & index_e <= corr_seg*corr_range & index_e ~= source);
                followers_corr_diag{i}(j) = length(idx_corr_diag);

                idx_corr_subdiag = find(index_e >= mod(1+corr_seg*corr_range, Ne) & ...
                        index_e <= mod((corr_seg+1)*corr_range-1, Ne)+1);
                followers_corr_subdiag{i}(j) = length(idx_corr_subdiag);
                
                if ~isempty(idx_corr_subdiag)
                        idx_corr_subdiag2 = find(index_e >= mod(1+(corr_seg+1)*corr_range, Ne) & ...
                                index_e <= mod((corr_seg+2)*corr_range-1, Ne)+1);
                        pre_neuron = idx_corr_subdiag;
                        flag = 0;
                        while flag == 0 && length(intersect(pre_neuron, idx_corr_subdiag2)) < length(idx_corr_subdiag2)
                                flag = 1;
                                for k = 1:length(idx_corr_subdiag2)
                                        if max(effw(index_e(idx_corr_subdiag2(k)), index_e)) == ...
                                                        max(effw(index_e(idx_corr_subdiag2(k)), index_e(pre_neuron))) && ...
                                                        sum(pre_neuron == idx_corr_subdiag2(k)) == 0
                                                pre_neuron = [pre_neuron; idx_corr_subdiag2(k)];
                                                flag = 0;
                                        end
                                end
                        end
                        idx_corr_subdiag2 = intersect(idx_corr_subdiag2, pre_neuron);
                else
                        idx_corr_subdiag2 = [];
                end
                followers_corr_subdiag2{i}(j) = length(idx_corr_subdiag2);
                
                idx_uncorr = setdiff(1:length(index_e), union(union(union(idx_corr_diag, idx_corr_subdiag), ...
                        idx_corr_subdiag2), find(index_e == source)))';
                followers_uncorr{i}(j) = length(idx_uncorr);
                
                avg_prob_corr_diag{i}(j) = mean(prob_follower_e(idx_corr_diag));
                avg_prob_corr_subdiag{i}(j) = mean(prob_follower_e(idx_corr_subdiag));
                avg_prob_corr_subdiag2{i}(j) = mean(prob_follower_e(idx_corr_subdiag2));
                avg_prob_uncorr{i}(j) = mean(prob_follower_e(idx_uncorr));
                all_prob_corr_diag{i} = [all_prob_corr_diag{i}; prob_follower_e(idx_corr_diag)];
                all_prob_corr_subdiag{i} = [all_prob_corr_subdiag{i}; prob_follower_e(idx_corr_subdiag)];
                all_prob_corr_subdiag2{i} = [all_prob_corr_subdiag2{i}; prob_follower_e(idx_corr_subdiag2)];
                all_prob_uncorr{i} = [all_prob_uncorr{i}; prob_follower_e(idx_uncorr)];
                avg_delay_corr_diag{i}(j) = mean(med_delay_e(idx_corr_diag));
                avg_delay_corr_subdiag{i}(j) = mean(med_delay_e(idx_corr_subdiag));
                avg_prob_corr_subdiag2{i}(j) = mean(med_delay_e(idx_corr_subdiag2));
                avg_delay_uncorr{i}(j) = mean(med_delay_e(idx_uncorr));
                all_delay_corr_diag{i} = [all_delay_corr_diag{i}; med_delay_e(idx_corr_diag)];
                all_delay_corr_subdiag{i} = [all_delay_corr_subdiag{i}; med_delay_e(idx_corr_subdiag)];
                all_delay_corr_subdiag2{i} = [all_delay_corr_subdiag2{i}; med_delay_e(idx_corr_subdiag2)];
                all_delay_uncorr{i} = [all_delay_uncorr{i}; med_delay_e(idx_uncorr)];
                cd('..');
        end
        cd('..')
end

save('followers_corr_uncorr4.mat', 'followers_corr_diag', 'followers_corr_subdiag', 'followers_corr_subdiag2', ...
        'followers_uncorr', 'avg_prob_corr_diag', 'avg_prob_corr_subdiag', 'avg_prob_corr_subdiag2', 'avg_prob_uncorr', ...
        'all_prob_corr_diag', 'all_prob_corr_subdiag', 'all_prob_corr_subdiag2', 'all_prob_uncorr', ...
        'avg_delay_corr_diag', 'avg_delay_corr_subdiag', 'avg_delay_corr_subdiag2', 'avg_delay_uncorr', ...
        'all_delay_corr_diag', 'all_delay_corr_subdiag', 'all_delay_corr_subdiag2', 'all_delay_uncorr');


