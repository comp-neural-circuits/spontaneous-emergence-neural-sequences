Ne = 1200;
Ni = 240;
sim_num = 10;

groups = 15;

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

cd('blank')
index_e_blank = 1:Ne;
index_i_blank = 1:Ni;
for i = 1:sim_num
       load(['sequence_blank_', num2str(i), '.mat'])
       index_e_blank = intersect(index_e_blank, follower_index_e);
       index_i_blank = intersect(index_i_blank, follower_index_i);
end
cd('..')

followers_e_blank = length(index_e_blank);
followers_i_blank = length(index_i_blank);

save('followers.mat', 'followers_e', 'followers_i', 'followers_e_blank', 'followers_i_blank');

                
                