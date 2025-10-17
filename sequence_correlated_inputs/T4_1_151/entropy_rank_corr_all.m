folders = dir();

id_delete = zeros(length(folders), 1);
        
for id = 1:length(id_delete)
        if ~isstrprop(folders(id).name(1), 'digit')        
                id_delete(id) = 1;
        end
end

folders(id_delete == 1) = [];
        
for id = 1:length(folders)
        cd(folders(id).name)
        sub_folders = dir();
        sub_folders = sub_folders(3:end);
        for sub_id = 1:length(sub_folders)
                cd(sub_folders(sub_id).name);
                disp(['Group ', num2str(id), ', sequence ', num2str(sub_id)])
                delay_and_jitter;
                calc_entropy;
                rank_correlation;
                cd('..');
        end
        cd('..');
end

