seq_dir = dir('T*');

for i = 1:length(seq_dir)
        cd(seq_dir(i).name)
        find_seq_all;
        delay_and_jitter_all;
        cd('..')
end