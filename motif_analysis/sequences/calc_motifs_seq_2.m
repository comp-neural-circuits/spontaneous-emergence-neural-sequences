
%seq_dir = '../../clean_sequence/';
current_folder = pwd;
net_dir = '../../basic_network/1/';
seq_dir = '../../sequence/';

folders = dir(seq_dir);

id_delete = zeros(length(folders), 1);
        
for id = 1:length(id_delete)
        if ~isstrprop(folders(id).name(1), 'digit')        
                id_delete(id) = 1;
        end
end

folders(id_delete == 1) = [];

ckpt = 80;
cutoff = 0.8;
V0 = -0.063;
th0 = 2.5;

EE_seg = h5read_sparsematrix([net_dir, '1/all_w.h5'], [num2str(ckpt), '/EE']);
EE = retrieve_sp(EE_seg);
ET = h5read([net_dir, '1/all_w.h5'], ['/', num2str(ckpt), '/ET']);
effw = EE./((ET - V0)*ones(1,size(EE, 2)))/1e3;
%{
n_seq = zeros(155, 4);
n_fol_1 = zeros(155, 1);
n_fol_2 = zeros(155, 1);
seq = 1;
        
for id = 1:length(folders)
        cd([seq_dir, folders(id).name])
        sub_folders = dir();
        sub_folders = sub_folders(3:end);
        for sub_id = 1:length(sub_folders)
                disp([id, sub_id])
                load([sub_folders(sub_id).name, '/sequence/delay_and_jitter.mat']);
                
                n_fol_1(seq) = length(index_e);  %followers
                followers = index_e(prob_follower_e > cutoff);
                n_fol_2(seq) = length(followers);  %followers with prob > cutoff
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
                seq = seq + 1;
        end
        cd(current_folder)
end

save(['motifs_seq_prob_', num2str(cutoff), '.mat'], 'n_fol_1', 'n_fol_2', 'n_seq');
%}

load(['motifs_seq_prob_', num2str(cutoff), '.mat'])
load('../basic_network/triplet_motifs_EE_id.mat')

n_net = squeeze(n(ckpt+1, 1, 1:4))';
prob_net = n_net/sum(n_net);
prob_seq = n_seq./sum(n_seq, 2);

weight = sum(n_seq, 2);
avg_prob_seq = nansum(prob_seq .* weight)/sum(weight);
std_prob_seq = sqrt(nansum(weight .* (prob_seq - avg_prob_seq).^2)/sum(weight)/...
        (sum(weight > 0)-1)*sum(weight > 0));
t_prob_seq = (avg_prob_seq - prob_net)./nanstd(prob_seq)*sqrt(sum(weight > 0));
p_prob_seq = 0.5 - abs(0.5 - tcdf(t_prob_seq, sum(weight > 0)));

prob_seq_norm = prob_seq ./ prob_net;
avg_prob_seq_norm = nansum(prob_seq_norm .* weight)/sum(weight);
std_prob_seq_norm = sqrt(nansum(weight .* (prob_seq_norm - avg_prob_seq_norm).^2)/...
        sum(weight)/(sum(weight > 0)-1)*sum(weight > 0));
t_prob_seq_norm = (avg_prob_seq_norm - 1)./nanstd(prob_seq_norm)*sqrt(sum(weight > 0));
p_prob_seq_norm = 0.5 - abs(0.5 - tcdf(t_prob_seq_norm, sum(weight > 0)));

color_motif{1} = [46, 49, 146]/255;
color_motif{2} = [0, 146, 69]/255;
color_motif{3} = [247, 147, 30]/255;
color_motif{4} = [117, 76, 36]/255;

figure
hold on
for i = 1:4
        scatter(i + 0.1*randn(sum(weight>0), 1), prob_seq_norm(weight>0, i), weight(weight>0), ...
                'MarkerEdgeColor', 'none', 'MarkerFaceColor', color_motif{i}, 'MarkerFaceAlpha', 0.6)
end
errorbar(avg_prob_seq_norm, std_prob_seq_norm, 'k', 'LineWidth', 2, 'LineStyle', 'none')
plot(avg_prob_seq_norm, 'ko', 'LineWidth', 2)
plot([0.5 4.5], [1 1], 'k--', 'LineWidth', 1)
text(1, 1.8, '*', 'FontSize', 20, 'HorizontalAlignment', 'center')
%text(2, 3.5, '*', 'FontSize', 20, 'HorizontalAlignment', 'center')
text(3, 1.5, '***', 'FontSize', 20, 'HorizontalAlignment', 'center')
text(4, 3.8, '***', 'FontSize', 20, 'HorizontalAlignment', 'center')
xlim([0.5 4.5])
ylim([0 4])
xticks(1:4)
xticklabels({})
set(gcf, 'color', [1 1 1])
ylabel('% of motif in sequences (normalized)')
yticks(0:4)
ax = gca;
a = get(ax,'Label');
set(gca,'Label', a,'fontsize',10)
set(gca, 'Box', 'off')
set(gca, 'TickDir', 'out')
savefig(gcf, 'figures/motifs_percentage_ratio.fig')
saveas(gcf, 'figures/motifs_percentage_ratio.pdf')

%load('../../feedback_inhibition/source/testing_inject_point_143.mat');
load('../../clean_sequence/testing_inject_point_143.mat');
figure
set(gcf, 'color', [1 1 1])

subplot(2, 2, 1)
plot(-100:0.1:300, rcd.Ve(1063, :)*1e3, 'LineWidth', 2)
xlim([0, 6])

subplot(2, 2, 2)
plot(-100:0.1:300, rcd.Ve(166, :)*1e3, 'LineWidth', 2)
xlim([6, 10])

subplot(2, 2, 3)
plot(-100:0.1:300, rcd.Ve(172, :)*1e3, 'LineWidth', 2)
xlim([7, 17])

subplot(2, 2, 4)
plot(-100:0.1:300, rcd.Ve(650, :)*1e3, 'LineWidth', 2)
xlim([5, 22])
savefig(gcf, 'figures/double_epsp.fig')
saveas(gcf, 'figures/double_epsp.pdf')


