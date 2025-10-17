file_name = 'test_net.jl';

template_dir = 'template';

groups = 19;

for i = 1:groups
        cd(num2str(i))
        folders = dir();
        folders = folders(3:end);
        for j = 1:length(folders)
                copyfile(['../', template_dir, '/', file_name], folders(j).name)
        end
        cd('..')
end

