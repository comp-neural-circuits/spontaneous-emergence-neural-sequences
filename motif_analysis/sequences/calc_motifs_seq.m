%seq_dir = '../../clean_sequence/';
net_dir = '../../basic_network/1/';
seq_dir = '../../feedback_inhibition/source/';

file = dir([seq_dir, '*.mat']);

ckpt = 80;
cutoff = 10;
V0 = -0.063;
th0 = 2.5;

n_seq = zeros(length(file), 4);
n_fol = zeros(length(file), 1);

for seq = 1:length(file)
        load([seq_dir, file(seq).name]);
        spikes = find(sum(rcd.xe, 2)>0);
        spike_time = zeros(length(spikes), 1);
        for i = 1:length(spikes)
                spike_time(i) = find(rcd.xe(spikes(i), :) == 1, 1);
        end
        M = [spikes, spike_time];
        sortM = sortrows(M, 2);
        n_fol(seq) = length(spikes);
        if cutoff == -1 || length(sortM) < cutoff
                followers = sortM(:, 1);
        else
                followers = sortM(1:cutoff, 1);
        end
        EE_seg = h5read_sparsematrix([net_dir, '1/all_w.h5'], [num2str(ckpt), '/EE']);
        EE = retrieve_sp(EE_seg);
        ET = h5read([net_dir, '1/all_w.h5'], ['/', num2str(ckpt), '/ET']);
        effw = EE./((ET - V0)*ones(1,size(EE, 2)))/1e3;
        backbone = (effw(followers, followers) > th0);
        
        for i = 1:length(followers)
                for j = 1:length(followers)
                        if backbone(j, i) == 1 && j ~= i
                                for k = 1:length(followers)
                                        if backbone(k, j) == 1 && k ~= i && k~= j
                                                if backbone(k, i) == 0 && backbone(i, j) == 0 && ...
                                                                backbone(i, k) == 0 && backbone(j, k) == 0
                                                        n_seq(seq, 1) = n_seq(seq, 1)+1;
                                                end
                                        end
                                                
                                        if backbone(k, i) == 1 && k ~= i && k > j
                                                if backbone(i, j) == 0 && backbone(i, k) == 0 && ...
                                                                        backbone(j, k) == 0 && backbone(k, j) == 0
                                                            n_seq(seq, 2) = n_seq(seq, 2)+1;
                                                end
                                        end
                                                
                                        if backbone(k, i) == 1 && backbone(k, j) == 1 && k~= i && k ~= j
                                                if backbone(i, j) == 0 && backbone(i, k) == 0 && backbone(j, k) == 0
                                                        n_seq(seq, 4) = n_seq(seq, 4)+1;
                                                end
                                        end
                                                
                                end
                        end
                                
                        if backbone(i, j) == 1 && j ~= i
                                for k = 1:length(followers)
                                        if backbone(i, k) == 1 && k ~= i && k > j
                                                if backbone(j, i) == 0 && backbone(k, i) == 0 && ...
                                                                        backbone(j, k) == 0 && backbone(k, j) == 0
                                                        n_seq(seq, 3) = n_seq(seq, 3)+1;
                                                end
                                        end
                                end
                        end
                end
        end
end

load('../basic_network/triplet_motifs_EE_id.mat')

n_net = squeeze(n(ckpt+1, 1, 1:4))';
%prob_net = n_net/sum(n_net);
%prob_seq = n_seq./sum(n_seq, 2);
%weights = sum(n_seq, 2);
%avg_prob_seq = sum(n_seq)/sum(sum(n_seq));
%std_prob_seq = sqrt(sum(weights .* (prob_seq - avg_prob_seq).^2)/sum(weights)*length(weights)/(length(weights)-1));
prob_net = n_net/1200;
prob_seq = n_seq./min(10, n_fol);

