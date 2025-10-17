load('followers.mat')
%{
if ~isfolder('delay_and_jitter')
        mkdir('delay_and_jitter')
end
%}
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

delay_all_e = [];
jitter_all_e = [];
prob_all_e = [];
delay_all_i = [];
jitter_all_i = [];
prob_all_i = [];
for i = 1:length(samples)
        load([num2str(sample_groups(i)), '/', num2str(samples(i)), '/sequence/delay_and_jitter.mat']);
        delay_all_e = [delay_all_e; med_delay_e];
        jitter_all_e = [jitter_all_e; jitter_e];
        prob_all_e = [prob_all_e; prob_follower_e];
        delay_all_i = [delay_all_i; med_delay_i];
        jitter_all_i = [jitter_all_i; jitter_i];
        prob_all_i = [prob_all_i; prob_follower_i]; 
end

save('delay_and_jitter_all.mat', 'delay_all_e', 'jitter_all_e', 'prob_all_e', 'delay_all_i', 'jitter_all_i', 'prob_all_i');
   