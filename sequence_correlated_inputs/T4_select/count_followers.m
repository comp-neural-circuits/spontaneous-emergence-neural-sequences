Ne = 1200;
Ni = 240;
sim_num = 10;

groups = 1;

followers_e = cell(groups, 1);
followers_i = cell(groups, 1);

for group = 1:groups
        cd(num2str(group))
        folders = dir();
        folders = folders(3:end);
        followers_e{group, 1} = zeros(length(folders), 1);
        followers_i{group, 1} = zeros(length(folders), 1);
        for site = 1:length(folders)
                cd(folders(site).name)
                index_e = 1:Ne;
                index_i = 1:Ni;
                for i = 1:sim_num
                        load(['sequence_', num2str(i), '.mat'])
                        index_e = intersect(index_e, follower_index_e);
                        index_i = intersect(index_i, follower_index_i);
                end
                followers_e{group, 1}(site) = length(index_e) - sum(index_e == str2double(folders(site).name));
                followers_i{group, 1}(site) = length(index_i);
                cd('..');
        end
        cd('..');
end

save('followers.mat', 'followers_e', 'followers_i');

                
                