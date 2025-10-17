folder_name = pwd;
pos = find(folder_name =='\', 1, 'last');
ckpt = [80:100, 110:10:160];

temp_dir = 'template';

for j = 1:length(ckpt)
        new_dir = [folder_name(pos+2:end), '_', num2str(ckpt(j))];
        mkdir(new_dir);
        copyfile(temp_dir, new_dir);
end


        