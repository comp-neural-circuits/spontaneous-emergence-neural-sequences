groups = 1;

for group = 1:groups
        cd(num2str(group))
        folders = dir();
        folders = folders(3:end);
        for site = 1:length(folders)
                disp(['Processing group ', num2str(group), ', neuron ', num2str(site)])
                cd(folders(site).name)
                find_seq;
                cd('..');
        end
        cd('..')
end