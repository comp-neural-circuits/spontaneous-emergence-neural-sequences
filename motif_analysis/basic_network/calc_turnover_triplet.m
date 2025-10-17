motif_dir = '';

resolution = [100, 200, 500, 1000, 1500, 2000];

ckpt_num = 21;     %total number of checkpoints

trials = 10;

motif_types = 5;

start_point = 8000;   %starting point of new net
early_interval_w = 100;  %time interval to save weights for the old net
ckpt_offset = 100;
interval_w = 100;

load([motif_dir, 'triplet_motifs_EE_id.mat']);

n0 = cell(trials, length(resolution));
n1 = cell(trials, length(resolution));
n2 = cell(trials, length(resolution));

for r = 1:length(resolution)
        all_ckpt = (1:round(resolution(r)/interval_w):min(ckpt_num, 1 + 10*round(resolution(r)/interval_w))) + ...
                round(start_point/early_interval_w);
        for trial = 1:trials
                disp(['Resolution: ', num2str(resolution(r)), 's, trial ' num2str(trial)])
                all_motif_id = motif_id{trial}(:, all_ckpt);
                n0{trial, r} = zeros(length(all_ckpt) - 1, motif_types);     %same
                n1{trial, r} = zeros(length(all_ckpt) - 1, motif_types);     %disappeared
                n2{trial, r} = zeros(length(all_ckpt) - 1, motif_types);     %emerged
                for motif = 1:motif_types
                        for j = 1:(length(all_ckpt) - 1)
                               if  isempty(all_motif_id{motif, j})
                                       n0{trial, r}(j, motif) = 0;
                                       n1{trial, r}(j, motif) = 0;
                                       n2{trial, r}(j, motif) = size(all_motif_id{motif, j+1}, 1);
                               else
                                       if isempty(all_motif_id{motif, j+1})
                                               n0{trial, r}(j, motif) = 0;
                                               n1{trial, r}(j, motif) = size(all_motif_id{motif, j}, 1);
                                               n2{trial, r}(j, motif) = 0;
                                       else
                                               n0{trial, r}(j, motif) = size(intersect(all_motif_id{motif, j}, all_motif_id{motif, j+1}, 'row'), 1);
                                               n1{trial, r}(j, motif) = size(all_motif_id{motif, j}, 1) - n0{trial, r}(j, motif);
                                               n2{trial, r}(j, motif) = size(all_motif_id{motif, j+1}, 1) - n0{trial, r}(j, motif);
                                       end
                               end
                        end
                end
        end
end

save('triplet_motif_turnover.mat', 'resolution', 'n0', 'n1', 'n2');

        
        
