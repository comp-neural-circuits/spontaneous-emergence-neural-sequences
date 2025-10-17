folders = dir();

id_delete = zeros(length(folders), 1);
        
for id = 1:length(id_delete)
        if ~isstrprop(folders(id).name(1), 'digit')        
                id_delete(id) = 1;
        elseif isfolder([folders(id).name, '/sequence'])
                id_delete(id) = 1;
        end
end

folders(id_delete == 1) = [];
        
for id = 1:length(folders)
        disp(['Processing sequence ', num2str(id)])
        cd(folders(id).name)
        if ~isfolder('sequence')
                delay_and_jitter;
        end
        cd('..');
end