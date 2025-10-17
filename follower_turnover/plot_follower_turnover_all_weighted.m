trials = 10;
source_neurons = cell(trials, 1);

source_neurons{1} = [1, 3, 5, 6, 9];
source_neurons{2} = [1, 3, 4, 5, 6];
source_neurons{3} = [1, 2, 3, 5, 6];
source_neurons{4} = [1, 2, 5, 6, 7];
source_neurons{5} = [1, 4, 6, 8, 9];
source_neurons{6} = [2, 3, 4, 5, 6];
source_neurons{7} = [1, 4, 5, 6, 7];
source_neurons{8} = [1, 3, 4, 5, 6];
source_neurons{9} = [1, 2, 4, 5, 6];
source_neurons{10} = [1, 2, 5, 6, 7];

net_1_folder = '../basic_network/2/';
net_2_folder = '../basic_network/3/';

if ~isfolder('figures')
        mkdir('figures')
end

ckpt_1 = [1:21, 31:10:71];
ckpt_2 = 1100;

ckpt_num1 = 71;
ckpt_num2 = 1001;

resolution1 = 100;
resolution2 = 1;
start_time = 8000;
offset1 = 100;
offset2 = 100;

time_1 = ckpt_1*resolution1 - offset1 + start_time;
time_2 = ckpt_2*resolution2 + ckpt_num1*resolution1 - offset1 - offset2 + start_time;
time_all = [time_1, time_2];

prob_bin = 0.1;                %the bin of the prob of followers
num_prob_bin = round(1/prob_bin);

max_delay = 130e-3;    %130ms, assumed maximal delay
delay_bin = 10e-3;         %10ms, the bin of delay
num_delay_bin = round(max_delay/delay_bin);

distance_bin = 1;                %the bin of the prob of followers
max_distance = 10;
num_distance_bin = round(max_distance/distance_bin);
Ne = 1200;
th0 = 2.5;

source_idx_all = cell2mat(source_neurons)';
trial_idx_all = ones(5, 1) * (1:10);
source_idx_all = source_idx_all(:);
trial_idx_all = trial_idx_all(:);
V0 = -0.063;

effw_1 = zeros(ckpt_num1, length(source_idx_all));
time_effw_1 = resolution1*(1:ckpt_num1) - offset1 + start_time;
effw_2 = zeros(round(ckpt_num2/resolution1), length(source_idx_all));
time_effw_2 = resolution2*(resolution1:resolution1:ckpt_num2) + ckpt_num1*resolution1 - offset1 -offset2 + start_time;

%{

num_follower_e_1 = zeros(length(ckpt_1), length(source_idx_all));
num_follower_i_1 = zeros(length(ckpt_1), length(source_idx_all));
num_follower_e_2 = zeros(length(ckpt_2), length(source_idx_all));
num_follower_i_2 = zeros(length(ckpt_2), length(source_idx_all));

index_e_all_1 = cell(length(ckpt_1), length(source_idx_all));
index_i_all_1 = cell(length(ckpt_1), length(source_idx_all));
index_e_all_2 = cell(length(ckpt_2), length(source_idx_all));
index_i_all_2 = cell(length(ckpt_2), length(source_idx_all));

prob_e_all_1 = cell(length(ckpt_1), length(source_idx_all));
prob_i_all_1 = cell(length(ckpt_1), length(source_idx_all));
prob_e_all_2 = cell(length(ckpt_2), length(source_idx_all));
prob_i_all_2 = cell(length(ckpt_2), length(source_idx_all));

delay_e_all_1 = cell(length(ckpt_1), length(source_idx_all));
delay_i_all_1 = cell(length(ckpt_1), length(source_idx_all));
delay_e_all_2 = cell(length(ckpt_2), length(source_idx_all));
delay_i_all_2 = cell(length(ckpt_2), length(source_idx_all));

distance_e_all_1 = cell(length(ckpt_1), length(source_idx_all));
distance_e_all_2 = cell(length(ckpt_2), length(source_idx_all));

for source = 1:length(source_idx_all)
        
        source_idx = source_idx_all(source);
        trial_idx = trial_idx_all(source);
        
        seq_dir = ['T', num2str(trial_idx), '_', num2str(source_idx), '/'];


        %seq_1_dir = ['T6_7_aeif_', num2str(trial_idx), '_', num2str(source_idx), '/'];
        %seq_2_dir = ['T6_8_aeif_', num2str(trial_idx), '_', num2str(source_idx), '/'];
        
        net_1_dir = [net_1_folder, num2str(trial_idx), '/'];
        net_2_dir = [net_2_folder, num2str(trial_idx), '/'];

        for i = 1:ckpt_num1
                EE_seg = h5read_sparsematrix([net_1_dir, 'all_w.h5'], [num2str(i), '/EE']);
                ET = h5read([net_1_dir, 'all_w.h5'], ['/', num2str(i), '/ET']);
                EE = retrieve_sp(EE_seg);
                effw = max(EE./(ET - V0))/1e3;
                effw_1(i, source) = effw(source_idx);
        end

        for i = 1:round(ckpt_num2/resolution1)
                EE_seg = h5read_sparsematrix([net_2_dir, 'all_w.h5'], [num2str(resolution1*i), '/EE']);
                ET = h5read([net_2_dir, 'all_w.h5'], ['/', num2str(resolution1*i), '/ET']);
                EE = retrieve_sp(EE_seg);
                effw = max(EE./(ET - V0))/1e3;
                effw_2(i, source) = effw(source_idx);
        end

        for i = 1:length(ckpt_1)
                load([seq_dir, num2str(trial_idx), '_', num2str(source_idx), '_', num2str(ckpt_1(i)+round(start_time/resolution1)-1), ...
                        '/sequence/delay_and_jitter.mat']);
                num_follower_e_1(i, source) = length(index_e) - 1;
                num_follower_i_1(i, source) = length(index_i);
                index_e_all_1{i, source} = setdiff(index_e, source_idx);
                index_i_all_1{i, source} = index_i;
                prob_e_all_1{i, source} = prob_follower_e(index_e ~= source_idx);
                prob_i_all_1{i, source} = prob_follower_i;
                delay_e_all_1{i, source} = med_delay_e(index_e ~= source_idx);
                delay_i_all_1{i, source} = med_delay_i;
                EE_seg = h5read_sparsematrix([net_1_dir, 'all_w.h5'], [num2str(ckpt_1(i)), '/EE']);
                ET = h5read([net_1_dir, 'all_w.h5'], ['/', num2str(ckpt_1(i)), '/ET']);
                EE = retrieve_sp(EE_seg);
                effw = EE./((ET - V0)*ones(1, Ne))/1e3;
                effw_backbone = digraph(effw' > th0);
                index_e_all_1{i, source} = setdiff(index_e, source_idx);
                distance_e_all_1{i, source} = zeros(size(index_e_all_1{i, source}));
                for j = 1:length(distance_e_all_1{i, source})
                       [path, distance] = shortestpath(effw_backbone, source_idx, index_e_all_1{i, source}(j));
                       distance_e_all_1{i, source}(j) = distance;
                end
                index_i_all_1{i, source} = index_i;
        end
        for i = 1:length(ckpt_2)
                load([seq_dir, num2str(trial_idx), '_', num2str(source_idx), '_', ...
                        num2str(round(ckpt_2(i)*resolution2 - offset2 + start_time)/resolution1 + ckpt_num1 - 1), ...
                        '/sequence/delay_and_jitter.mat']);
                num_follower_e_2(i, source) = length(index_e) - 1;
                num_follower_i_2(i, source) = length(index_i);
                index_e_all_2{i, source} = setdiff(index_e, source_idx);
                index_i_all_2{i, source} = index_i;
                prob_e_all_2{i, source} = prob_follower_e(index_e ~= source_idx);
                prob_i_all_2{i, source} = prob_follower_i;
                delay_e_all_2{i, source} = med_delay_e;
                delay_i_all_2{i, source} = med_delay_i;
                EE_seg = h5read_sparsematrix([net_2_dir, 'all_w.h5'], [num2str(ckpt_2(i)), '/EE']);
                ET = h5read([net_2_dir, 'all_w.h5'], ['/', num2str(ckpt_2(i)), '/ET']);
                EE = retrieve_sp(EE_seg);
                effw = EE./((ET - V0)*ones(1, Ne))/1e3;
                effw_backbone = digraph(effw' > th0);
                index_e_all_2{i, source} = setdiff(index_e, source_idx);
                distance_e_all_2{i, source} = zeros(size(index_e_all_2{i, source}));
                for j = 1:length(distance_e_all_2{i, source})
                        [path, distance] = shortestpath(effw_backbone, source_idx, index_e_all_2{i, source}(j));
                        distance_e_all_2{i, source}(j) = distance;
                end
                index_i_all_2{i, source} = index_i;
        end
end

save('follower_turnover_data.mat', 'effw_1', 'time_effw_1', 'effw_2', 'time_effw_2', ...
        'num_follower_e_1', 'num_follower_i_1', 'num_follower_e_2', 'num_follower_i_2', ...
        'index_e_all_1', 'index_i_all_1', 'index_e_all_2', 'index_i_all_2', ...
        'prob_e_all_1', 'prob_i_all_1', 'prob_e_all_2', 'prob_i_all_2', ...
        'delay_e_all_1', 'delay_i_all_1', 'delay_e_all_2', 'delay_i_all_2', 'distance_e_all_1', 'distance_e_all_2', ...
        'source_idx_all', 'trial_idx_all');

%}

load('follower_turnover_data.mat');

index_e_all_3 = index_e_all_1(1:21, :);   %from 8000 to 10000
prob_e_all_3 = prob_e_all_1(1:21, :);
delay_e_all_3 = delay_e_all_1(1:21, :);
distance_e_all_3 = distance_e_all_1(1:21, :);
time_3 = time_all(1:21);
ratio1 = zeros(size(index_e_all_3, 1)-1, length(source_idx_all));   %100s interval, old followers
ratio2 = zeros(size(index_e_all_3, 1)-1, length(source_idx_all));   %100s interval, new followers
ratio3 = ones(size(index_e_all_3, 1), length(source_idx_all));        %compare to the first frame
w_ratio3 = ones(1, length(source_idx_all));
for source = 1:length(source_idx_all)
        w_ratio3(source) = length(index_e_all_3{1, source});
        for i = 1:(size(index_e_all_3, 1)-1)
                ratio1(i, source) = length(intersect(index_e_all_3{i, source}, index_e_all_3{i+1, source}))/length(index_e_all_3{i, source});
                ratio2(i, source) = length(intersect(index_e_all_3{i, source}, index_e_all_3{i+1, source}))/length(index_e_all_3{i+1, source});
                ratio3(i+1, source) = length(intersect(index_e_all_3{1, source}, index_e_all_3{i+1, source}))/length(index_e_all_3{1, source});
        end
end

avg_ratio3 = sum(ratio3' .* w_ratio3')'/sum(w_ratio3);
std_ratio3 = std_weight(ratio3', w_ratio3')';
figure
fill([time_3, flip(time_3)], [avg_ratio3 - std_ratio3; flip(avg_ratio3 + std_ratio3)], ...
        [0.7 0.7 1], 'LineStyle', 'none');
hold on
plot(time_3, avg_ratio3, 'b', 'LineWidth', 2);
f = fit(time_3'-min(time_3), avg_ratio3, '(1-c)*exp(-x/b)+c', 'StartPoint', [500, 0.2]);
text(9200, 0.8, ['\tau_{seq} = ', num2str(round(f.b)), 's'], 'color', 'b', 'FontSize', 15)
text(9200, 0.7, ['baseline = ', num2str(round(f.c, 2))], 'color', 'b', 'FontSize', 15)
xlabel('Time (s)')
ylabel('Ratio of remaining followers')
set(gcf, 'color', [1 1 1])
ax = gca;
a = get(ax,'Label');
set(gca,'Label', a,'fontsize',15)
set(gca, 'Box', 'off')
ylim([0 1])
savefig(gcf, 'figures/turnover_rate.fig');

figure
plot([time_effw_1, time_effw_2], [effw_1; effw_2], 'k')
hold on
plot([time_effw_1, time_effw_2], mean([effw_1; effw_2], 2), 'k', 'LineWidth', 3)
xlim([8000 16000])
xlabel('Time (s)')
ylabel('Weight-threshold ratio')
set(gcf, 'color', [1 1 1])
ax = gca;
a = get(ax,'Label');
set(gca,'Label', a,'fontsize',15)
set(gca, 'Box', 'off')
savefig(gcf, 'figures/effw.fig');

figure
plot(time_all, mean([num_follower_e_1; num_follower_e_2], 2), 'b', 'LineWidth', 3);
hold on
plot(time_all, mean([num_follower_i_1; num_follower_i_2], 2), 'r', 'LineWidth', 3);
plot(time_all, [num_follower_e_1; num_follower_e_2], 'b');
plot(time_all, [num_follower_i_1; num_follower_i_2], 'r');
xlabel('Time (s)')
ylabel('Number of followers')
legend('Excitatory', 'Inhibitory', 'Location', 'northwest')
set(gcf, 'color', [1 1 1])
ax = gca;
a = get(ax,'Label');
set(gca,'Label', a,'fontsize',15)
set(gca, 'Box', 'off')
savefig(gcf, 'figures/num_followers.fig');


index_e_all_4 = index_e_all_1([1:10:21, 22:end], :);   %from 8000 to 16000
prob_e_all_4 = prob_e_all_1([1:10:21, 22:end], :);
delay_e_all_4 = delay_e_all_1([1:10:21, 22:end], :);
distance_e_all_4 = distance_e_all_1([1:10:21, 22:end], :);
time_4 = time_all([1:10:21, 22:end]);

ratio_remain_e = zeros(num_prob_bin, length(source_idx_all));
w_remain_e = zeros(num_prob_bin, length(source_idx_all));
for source = 1:length(source_idx_all)
         for i = 1:num_prob_bin
                 w_k = zeros(size(index_e_all_4, 1)-1, 1);
                 ratio_k = zeros(size(index_e_all_4, 1)-1, 1);
                 for k = 1:size(index_e_all_4, 1)-1
                         w_k(k) = length(index_e_all_4{k, source}(prob_e_all_4{k, source} > (i-1)*prob_bin & ...
                                 prob_e_all_4{k, source} <= i*prob_bin));
                         if w_k(k) ~= 0
                                 ratio_k(k) = length(intersect(index_e_all_4{k, source}(prob_e_all_4{k, source} > (i-1)*prob_bin & ...
                                         prob_e_all_4{k, source} <= i*prob_bin), index_e_all_4{k+1, source})) / w_k(k);
                         end     
                 end
                 ratio_remain_e(i, source) = sum(ratio_k .* w_k)/sum(w_k);
                 w_remain_e(i, source) = sum(w_k);
         end
end

avg_ratio_remain_e = zeros(num_prob_bin, 1);
min_ratio_remain_e = zeros(num_prob_bin, 1);
max_ratio_remain_e = zeros(num_prob_bin, 1);
sem_ratio_remain_e = zeros(num_prob_bin, 1);
std_ratio_remain_e = zeros(num_prob_bin, 1);
for i = 1:num_prob_bin
        valid_idx = ~isnan(ratio_remain_e(i, :));
        avg_ratio_remain_e(i) = sum(ratio_remain_e(i, valid_idx) .* w_remain_e(i, valid_idx))/sum(w_remain_e(i, valid_idx));
        min_ratio_remain_e(i) = min(ratio_remain_e(i, valid_idx));
        max_ratio_remain_e(i) = max(ratio_remain_e(i, valid_idx));
        sem_ratio_remain_e(i) = sem_weight(ratio_remain_e(i, valid_idx)', w_remain_e(i, valid_idx)');
        std_ratio_remain_e(i) = std_weight(ratio_remain_e(i, valid_idx)', w_remain_e(i, valid_idx)');
end

figure
fill([(prob_bin/2):prob_bin:(1-prob_bin/2), flip((prob_bin/2):prob_bin:(1-prob_bin/2))], ...
        [(avg_ratio_remain_e - std_ratio_remain_e)'...
        , flip(avg_ratio_remain_e + std_ratio_remain_e)'], [101 200 208]/255, 'LineStyle', 'none');
hold on
fill([(prob_bin/2):prob_bin:(1-prob_bin/2), flip((prob_bin/2):prob_bin:(1-prob_bin/2))], ...
        [(avg_ratio_remain_e - sem_ratio_remain_e)'...
        , flip(avg_ratio_remain_e + sem_ratio_remain_e)'], [0.7 0.7 1], 'LineStyle', 'none');
plot((prob_bin/2):prob_bin:(1-prob_bin/2), avg_ratio_remain_e, 'b', 'LineWidth', 2);
legend('Std', 's.e.m', 'Average', 'Location', 'northwest')
xlabel('Probability in the initial sequence')
ylabel('Remaining ratio')
set(gcf, 'color', [1 1 1])
ax = gca;
a = get(ax,'Label');
set(gca,'Label', a,'fontsize',15)
set(gca, 'box', 'off')
ylim([0, ax.YLim(2)])
savefig(gcf, 'figures/turnover_vs_prob.fig');

ratio_remain_e_delay = zeros(num_delay_bin+1, length(source_idx_all));
w_remain_e_delay = zeros(num_delay_bin+1, length(source_idx_all));
for source = 1:length(source_idx_all)
         for i = 1:num_delay_bin
                w_k = zeros(11, 1);
                ratio_k = zeros(11, 1);
                for k = 1:size(index_e_all_4, 1)-1
                        w_k(k) = length(index_e_all_4{k, source}(delay_e_all_4{k, source} > ...
                                (i-1)*delay_bin & delay_e_all_4{k, source} <= i*delay_bin));
                        if w_k(k) ~= 0
                                ratio_k(k) = length(intersect(index_e_all_4{k, source}(delay_e_all_4{k, source} > ...
                                        (i-1)*delay_bin & delay_e_all_4{k, source} <= i*delay_bin), index_e_all_4{k+1, source})) / w_k(k);
                        end
                end
                ratio_remain_e_delay(i, source) = sum(ratio_k .* w_k)/sum(w_k);
                w_remain_e_delay(i, source) = sum(w_k);
         end
         w_k = zeros(size(index_e_all_4, 1)-1, 1);
         ratio_k = zeros(size(index_e_all_4, 1)-1, 1);
         for k = 1:size(index_e_all_4, 1)-1
                 w_k(k) = length(index_e_all_4{k, source}(delay_e_all_4{k, source} > max_delay));
                 if w_k(k) ~= 0
                         ratio_k(k) = length(intersect(index_e_all_4{k, source}...
                                 (delay_e_all_4{k, source} > max_delay), index_e_all_4{k+1, source}))/w_k(k);
                 end
         end
         ratio_remain_e_delay(num_delay_bin + 1, source) = sum(ratio_k .* w_k)/sum(w_k);
         w_remain_e_delay(num_delay_bin + 1, source) = sum(w_k);
end

avg_ratio_remain_e_delay = zeros(num_delay_bin + 1, 1);
min_ratio_remain_e_delay = zeros(num_delay_bin + 1, 1);
max_ratio_remain_e_delay = zeros(num_delay_bin + 1, 1);
sem_ratio_remain_e_delay = zeros(num_delay_bin + 1, 1);
std_ratio_remain_e_delay = zeros(num_delay_bin + 1, 1);
for i = 1:(num_delay_bin + 1)
        valid_idx = ~isnan(ratio_remain_e_delay(i, :));
        if sum(valid_idx) == 0
                avg_ratio_remain_e_delay(i) = NaN;
                min_ratio_remain_e_delay(i) = NaN;
                max_ratio_remain_e_delay(i) = NaN;
                sem_ratio_remain_e_delay(i) = NaN;
                std_ratio_remain_e_delay(i) = NaN;
        else
                avg_ratio_remain_e_delay(i) = sum(ratio_remain_e_delay(i, valid_idx) .* w_remain_e_delay(i, valid_idx))/...
                        sum(w_remain_e_delay(i, valid_idx));
                min_ratio_remain_e_delay(i) = min(ratio_remain_e_delay(i, valid_idx));
                max_ratio_remain_e_delay(i) = max(ratio_remain_e_delay(i, valid_idx));
                sem_ratio_remain_e_delay(i) = sem_weight(ratio_remain_e_delay(i, valid_idx)', w_remain_e_delay(i, valid_idx)');
                std_ratio_remain_e_delay(i) = std_weight(ratio_remain_e_delay(i, valid_idx)', w_remain_e_delay(i, valid_idx)');
        end
end

figure
fill([delay_bin/2:delay_bin:(max_delay+delay_bin/2), flip(delay_bin/2:delay_bin:(max_delay+delay_bin/2))]*1e3, ...
        [(avg_ratio_remain_e_delay - std_ratio_remain_e_delay)' ...
        , flip(avg_ratio_remain_e_delay + std_ratio_remain_e_delay)'], [101 200 208]/255, 'LineStyle', 'none');
hold on
fill([delay_bin/2:delay_bin:(max_delay+delay_bin/2), flip(delay_bin/2:delay_bin:(max_delay+delay_bin/2))]*1e3, ...
        [(avg_ratio_remain_e_delay - sem_ratio_remain_e_delay)' ...
        , flip(avg_ratio_remain_e_delay + sem_ratio_remain_e_delay)'], [0.7 0.7 1], 'LineStyle', 'none');
plot((delay_bin/2:delay_bin:(max_delay+delay_bin/2))*1e3, avg_ratio_remain_e_delay, 'b', 'LineWidth', 2);
legend('Std', 's.e.m', 'Average', 'Location', 'northeast')
xlabel('Median delay (ms)')
ylabel('Remaining ratio')
set(gcf, 'color', [1 1 1])
ax = gca;
a = get(ax,'Label');
set(gca,'Label', a,'fontsize',15)
set(gca, 'box', 'off')
ylim([0, ax.YLim(2)])
savefig(gcf, 'figures/turnover_vs_delay.fig');

ratio_remain_e_distance = zeros(num_distance_bin+1, length(source_idx_all));
w_remain_e_distance = zeros(num_distance_bin+1, length(source_idx_all));
for source = 1:length(source_idx_all)
         for i = 1:num_distance_bin
                w_k = zeros(size(index_e_all_4, 1)-1, 1);
                ratio_k = zeros(size(index_e_all_4, 1)-1, 1);
                for k = 1:size(index_e_all_4, 1)-1
                        w_k(k) = length(index_e_all_4{k, source}(distance_e_all_4{k, source} > ...
                                (i-1)*distance_bin & distance_e_all_4{k, source} <= i*distance_bin));
                        if w_k(k) ~= 0
                                ratio_k(k) = length(intersect(index_e_all_4{k, source}(distance_e_all_4{k, source} > ...
                                        (i-1)*distance_bin & distance_e_all_4{k, source} <= i*distance_bin), index_e_all_4{k+1, source})) / w_k(k);
                        end
                end
                ratio_remain_e_distance(i, source) = sum(ratio_k .* w_k)/sum(w_k);
                w_remain_e_distance(i, source) = sum(w_k);
         end
         w_k = zeros(size(index_e_all_4, 1)-1, 1);
         ratio_k = zeros(size(index_e_all_4, 1)-1, 1);
         for k = 1:size(index_e_all_4, 1)-1
                 w_k(k) = length(index_e_all_4{k, source}(distance_e_all_4{k, source} > max_distance));
                 if w_k(k) ~= 0
                         ratio_k(k) = length(intersect(index_e_all_4{k, source}...
                                 (distance_e_all_4{k, source} > max_distance), index_e_all_4{k+1, source}))/w_k(k);
                 end
         end
         ratio_remain_e_distance(num_distance_bin + 1, source) = sum(ratio_k .* w_k)/sum(w_k);
         w_remain_e_distance(num_distance_bin + 1, source) = sum(w_k);
end

avg_ratio_remain_e_distance = zeros(num_distance_bin+1, 1);
min_ratio_remain_e_distance = zeros(num_distance_bin+1, 1);
max_ratio_remain_e_distance = zeros(num_distance_bin+1, 1);
sem_ratio_remain_e_distance = zeros(num_distance_bin+1, 1);
std_ratio_remain_e_distance = zeros(num_distance_bin+1, 1);
for i = 1:(num_distance_bin+1)
        valid_idx = ~isnan(ratio_remain_e_distance(i, :));
        if sum(~isnan(ratio_remain_e_distance(i, :))) == 0
                avg_ratio_remain_e_distance(i) = NaN;
                min_ratio_remain_e_distance(i) = NaN;
                max_ratio_remain_e_distance(i) = NaN;
                sem_ratio_remain_e_distance(i) = NaN;
                std_ratio_remain_e_distance(i) = NaN;
        else
                avg_ratio_remain_e_distance(i) = sum(ratio_remain_e_distance(i, valid_idx) .* w_remain_e_distance(i, valid_idx))/...
                        sum(w_remain_e_distance(i, valid_idx));
                min_ratio_remain_e_distance(i) = min(ratio_remain_e_distance(i, valid_idx));
                max_ratio_remain_e_distance(i) = max(ratio_remain_e_distance(i, valid_idx));
                sem_ratio_remain_e_distance(i) = sem_weight(ratio_remain_e_distance(i, valid_idx)', w_remain_e_distance(i, valid_idx)');
                std_ratio_remain_e_distance(i) = std_weight(ratio_remain_e_distance(i, valid_idx)', w_remain_e_distance(i, valid_idx)');
        end
end

figure
fill([distance_bin:distance_bin:(max_distance+distance_bin), ...
        flip(distance_bin:distance_bin:(max_distance+distance_bin))], ...
        [(avg_ratio_remain_e_distance - std_ratio_remain_e_distance)'...
        , flip(avg_ratio_remain_e_distance + std_ratio_remain_e_distance)'], [101 200 208]/255, 'LineStyle', 'none');
hold on
fill([distance_bin:distance_bin:(max_distance+distance_bin), ...
        flip(distance_bin:distance_bin:(max_distance+distance_bin))], ...
        [(avg_ratio_remain_e_distance - sem_ratio_remain_e_distance)'...
        , flip(avg_ratio_remain_e_distance + sem_ratio_remain_e_distance)'], [0.7 0.7 1], 'LineStyle', 'none');
plot(distance_bin:distance_bin:(max_distance+distance_bin), avg_ratio_remain_e_distance, 'b', 'LineWidth', 2);
legend('Std', 's.e.m.', 'Average', 'Location', 'northeast')
xlabel('Distance in the initial connections')
ylabel('Remaining ratio')
set(gcf, 'color', [1 1 1])
ax = gca;
a = get(ax,'Label');
set(gca,'Label', a,'fontsize',15)
set(gca, 'box', 'off')
ylim([0, ax.YLim(2)])
savefig(gcf, 'figures/turnover_vs_distance.fig');

new_follower_e = zeros(size(index_e_all_3, 1) - 1, size(index_e_all_3, 2));
disappear_e = zeros(size(index_e_all_3, 1) - 1, size(index_e_all_3, 2));
for i = 1:(size(delay_e_all_3, 1) - 1)
        for j = 1:size(delay_e_all_3, 2)
                new_follower_e(i, j) = length(index_e_all_3{i+1, j}) - ...
                        length(intersect(index_e_all_3{i+1, j}, index_e_all_3{i, j}));
                disappear_e(i, j) = length(index_e_all_3{i, j}) - ...
                        length(intersect(index_e_all_3{i+1, j}, index_e_all_3{i, j}));
        end
end

figure
set(gcf, 'color', [1 1 1])
plot((time_3(1:end-1)+time_3(2:end))/2, mean(new_follower_e, 2), '-o', 'color', [247, 147, 30]/255, 'LineWidth', 2)
hold on
plot((time_3(1:end-1)+time_3(2:end))/2, mean(disappear_e, 2), '--o', 'color', 'k', 'LineWidth', 2)
xticks(8000:500:10000)
xticklabels({'8000', '8500', '9000', '9500', '10000'});
xlabel('Time (s)')
ylabel('Number of followers')
legend('New followers', 'Disappeared followers')
ax = gca;
a = get(ax,'Label');
set(gca,'Label', a,'fontsize',15)
set(gca, 'box', 'off')
savefig(gcf, 'figures/new_disappeared_followers.fig');

remain_e_prob = cell((size(index_e_all_3, 1) - 1)*size(index_e_all_3, 2), 1);
new_follower_e_prob = cell((size(index_e_all_3, 1) - 1)*size(index_e_all_3, 2), 1);
new_idx = 1;
remain_idx = 1;
for i = 1:(size(delay_e_all_3, 1) - 1)
        for j = 1:size(delay_e_all_3, 2)
                remain_e_idx = intersect(index_e_all_3{i+1, j}, index_e_all_3{i, j});
                remain_e_array_idx = zeros(size(remain_e_idx));
                for k = 1:length(remain_e_idx)
                        remain_e_array_idx(k) = find(remain_e_idx(k) == index_e_all_3{i+1, j});
                end
                new_follower_e_array_idx = setdiff(1:length(index_e_all_3{i+1, j}), remain_e_array_idx);
                if ~isempty(remain_e_array_idx)
                        remain_e_prob{remain_idx} = prob_e_all_3{i+1, j}(remain_e_array_idx);
                        remain_idx = remain_idx +1;
                end
                if ~isempty(new_follower_e_array_idx)
                        new_follower_e_prob{new_idx} = prob_e_all_3{i+1, j}(new_follower_e_array_idx);
                        new_idx = new_idx +1;
                end
        end
end
remain_e_prob_all = cell2mat(remain_e_prob(1:(remain_idx - 1)));
new_follower_e_prob_all = cell2mat(new_follower_e_prob(1:(new_idx - 1)));

figure
set(gcf, 'color', [1 1 1])
histogram(remain_e_prob_all, 'LineStyle', 'none', 'Normalization', 'pdf');
hold on
histogram(new_follower_e_prob_all, 'LineStyle', 'none', 'Normalization', 'pdf');
xlabel('Reliability of followers')
ylabel('Probability density')
legend('Remaining', 'New')
ax = gca;
a = get(ax,'Label');
set(gca,'Label', a,'fontsize',15)
set(gca, 'box', 'off')
savefig(gcf, 'figures/prob_new_disappeared_followers.fig');

remain_e_delay = cell((size(index_e_all_3, 1) - 1)*size(index_e_all_3, 2), 1);
new_follower_e_delay = cell((size(index_e_all_3, 1) - 1)*size(index_e_all_3, 2), 1);
new_idx = 1;
remain_idx = 1;
for i = 1:(size(delay_e_all_3, 1) - 1)
        for j = 1:size(delay_e_all_3, 2)
                remain_e_idx = intersect(index_e_all_3{i+1, j}, index_e_all_3{i, j});
                remain_e_array_idx = zeros(size(remain_e_idx));
                for k = 1:length(remain_e_idx)
                        remain_e_array_idx(k) = find(remain_e_idx(k) == index_e_all_3{i+1, j});
                end
                new_follower_e_array_idx = setdiff(1:length(index_e_all_3{i+1, j}), remain_e_array_idx);
                if ~isempty(remain_e_array_idx)
                        remain_e_delay{remain_idx} = delay_e_all_3{i+1, j}(remain_e_array_idx);
                        remain_idx = remain_idx +1;
                end
                if ~isempty(new_follower_e_array_idx)
                        new_follower_e_delay{new_idx} = delay_e_all_3{i+1, j}(new_follower_e_array_idx);
                        new_idx = new_idx +1;
                end
        end
end
remain_e_delay_all = cell2mat(remain_e_delay(1:(remain_idx - 1)));
new_follower_e_delay_all = cell2mat(new_follower_e_delay(1:(new_idx - 1)));

remain_e_distance = cell((size(index_e_all_3, 1) - 1)*size(index_e_all_3, 2), 1);
new_follower_e_distance = cell((size(index_e_all_3, 1) - 1)*size(index_e_all_3, 2), 1);
new_idx = 1;
remain_idx = 1;
for i = 1:(size(distance_e_all_3, 1) - 1)
        for j = 1:size(distance_e_all_3, 2)
                remain_e_idx = intersect(index_e_all_3{i+1, j}, index_e_all_3{i, j});
                remain_e_array_idx = zeros(size(remain_e_idx));
                for k = 1:length(remain_e_idx)
                        remain_e_array_idx(k) = find(remain_e_idx(k) == index_e_all_3{i+1, j});
                end
                new_follower_e_array_idx = setdiff(1:length(index_e_all_3{i+1, j}), remain_e_array_idx);
                if ~isempty(remain_e_array_idx)
                        remain_e_distance{remain_idx} = distance_e_all_3{i+1, j}(remain_e_array_idx);
                        remain_idx = remain_idx +1;
                end
                if ~isempty(new_follower_e_array_idx)
                        new_follower_e_distance{new_idx} = distance_e_all_3{i+1, j}(new_follower_e_array_idx);
                        new_idx = new_idx +1;
                end
        end
end
remain_e_distance_all = cell2mat(remain_e_distance(1:(remain_idx - 1)));
new_follower_e_distance_all = cell2mat(new_follower_e_distance(1:(new_idx - 1)));

all_followers_traces = cell(1000, 1);
all_followers_idx = zeros(1000, 1);
k = 1;
for i = 1:size(index_e_all_3, 2)
        index_non_empty = ones(size(index_e_all_3, 1), 1);
        for j = 1:length(index_non_empty)
                if isempty(index_e_all_3{j, i})
                       index_non_empty(j) = 0;
                end
        end
        all_followers = unique(cell2mat(index_e_all_3(index_non_empty == 1, i)));
        for j = 1:length(all_followers)
                all_followers_traces{k} = struct();
                all_followers_traces{k}.follower_idx = Ne*(i-1) + all_followers(j);
                all_followers_traces{k}.time = [];
                all_followers_traces{k}.prob = [];
                all_followers_idx(k) = Ne*(i-1) + all_followers(j);
                k = k + 1;
        end
end
        
all_followers_traces = all_followers_traces(1:k-1);

for i = 1:size(index_e_all_1, 2)
        for t = 1:size(index_e_all_1, 1)
                followers_idx = index_e_all_1{t, i};
                for j = 1:length(followers_idx)
                        idx = find(followers_idx(j) + Ne*(i-1) == all_followers_idx);
                        if ~isempty(idx)
                                all_followers_traces{idx}.time = [all_followers_traces{idx}.time; time_1(t)];
                                all_followers_traces{idx}.prob = [all_followers_traces{idx}.prob; ...
                                        prob_e_all_1{t, i}(all_followers_idx(idx) -  Ne*(i-1) == index_e_all_1{t, i})];
                        end
                end
        end
end

prob_all = cell(length(all_followers_idx), 1);
duration_all_new = zeros(length(all_followers_idx), 1);
rise_time_all_new = zeros(length(all_followers_idx), 1);
max_prob = zeros(length(all_followers_idx), 1);
init_prob = zeros(length(all_followers_idx), 1);

k = 1;
for i = 1:length(all_followers_idx)
        if min(all_followers_traces{i}.time) > start_time && max(all_followers_traces{i}.time) < time_3(end) %newly emerged followers
                %for j = 1:length(all_followers_traces{i}.time)
                    %    prob_all{j} = [prob_all{j}; all_followers_traces{i}.prob(j)];
                %end
                if length(all_followers_traces{i}.time) == round((all_followers_traces{i}.time(end) - ...
                                all_followers_traces{i}.time(1))/100) + 1
                        prob_all{k} = all_followers_traces{i}.prob;
                        duration_all_new(k) = length(all_followers_traces{i}.time);
                        rise_time_all_new(k) = find(all_followers_traces{i}.prob == max(all_followers_traces{i}.prob), 1);
                        max_prob(k) = max(all_followers_traces{i}.prob);
                        init_prob(k) = all_followers_traces{i}.prob(1);
                        k = k+1;
                end
        end
end
max_prob = max_prob(1:k-1);
init_prob = init_prob(1:k-1);
duration_all_new = duration_all_new(1:k-1);
rise_time_all_new = rise_time_all_new(1:k-1);
prob_all = prob_all(1:k-1);

prob_avg_traces = cell(max(duration_all_new), 1);
for i = 1:max(duration_all_new)
        prob_traces = cell2mat(prob_all(duration_all_new == i)');
        prob_avg_traces{i} = mean(prob_traces, 2);
end

figure
set(gcf, 'color', [1 1 1])
for i = 1:length(prob_avg_traces)
        if mod(i, 3) == 0
                plot((1:i)*100, prob_avg_traces{i}, '-o', 'LineWidth', 2)
                        %'color', [1 1 1]*0.8*(length(prob_avg_traces) - i)/(length(prob_avg_traces) - 1))
                hold on
        end
end
xlabel('Interval (s)')
ylabel('Average responding probability')
legend('300s', '600s', '900s', '1200s', '1500s', '1800s', 'Location', 'northwest')
ax = gca;
a = get(ax,'Label');
set(gca,'Label', a,'fontsize',15)
set(gca, 'box', 'off')
savefig(gcf, 'figures/follower_prob_traces.fig');

avg_max_prob = zeros(max(duration_all_new), 1);
avg_init_prob = zeros(max(duration_all_new), 1);
std_max_prob = zeros(max(duration_all_new), 1);
std_init_prob = zeros(max(duration_all_new), 1);
sem_max_prob = zeros(max(duration_all_new), 1);
sem_init_prob = zeros(max(duration_all_new), 1);

for i = 1:max(duration_all_new)
        avg_max_prob(i) = mean(max_prob(duration_all_new == i));
        std_max_prob(i) = std(max_prob(duration_all_new == i));
        sem_max_prob(i) = std_max_prob(i)/sqrt(sum(duration_all_new == i));
        avg_init_prob(i) = mean(init_prob(duration_all_new == i));
        std_init_prob(i) = std(init_prob(duration_all_new == i));
        sem_init_prob(i) = std_init_prob(i)/sqrt(sum(duration_all_new == i));
end

figure
fill([1:max(duration_all_new), flip(1:max(duration_all_new))]*resolution1, ...
        [(avg_max_prob - std_max_prob)', flip(avg_max_prob + std_max_prob)'], [101 200 208]/255, 'LineStyle', 'none');
hold on
fill([1:max(duration_all_new), flip(1:max(duration_all_new))]*resolution1, ...
        [(avg_max_prob - sem_max_prob)', flip(avg_max_prob + sem_max_prob)'], [0.7 0.7 1], 'LineStyle', 'none');
plot((1:max(duration_all_new))*resolution1, avg_max_prob, 'b', 'LineWidth', 2);
legend('Std', 's.e.m.', 'Average', 'Location', 'southeast')
xlabel('Interval (s)')
ylabel('Maximal responding probability')
set(gcf, 'color', [1 1 1])
ax = gca;
a = get(ax,'Label');
set(gca,'Label', a,'fontsize',15)
set(gca, 'box', 'off')
ylim([0, ax.YLim(2)])
savefig(gcf, 'figures/follower_duration_max_prob.fig');

figure
fill([1:max(duration_all_new), flip(1:max(duration_all_new))]*resolution1, ...
        [(avg_init_prob - std_init_prob)', flip(avg_init_prob + std_init_prob)'], [101 200 208]/255, 'LineStyle', 'none');
hold on
fill([1:max(duration_all_new), flip(1:max(duration_all_new))]*resolution1, ...
        [(avg_init_prob - sem_init_prob)', flip(avg_init_prob + sem_init_prob)'], [0.7 0.7 1], 'LineStyle', 'none');
plot((1:max(duration_all_new))*resolution1, avg_init_prob, 'b', 'LineWidth', 2);
legend('Std', 's.e.m.', 'Average', 'Location', 'southeast')
xlabel('Interval (s)')
ylabel('Maximal responding probability')
set(gcf, 'color', [1 1 1])
ax = gca;
a = get(ax,'Label');
set(gca,'Label', a,'fontsize',15)
set(gca, 'box', 'off')
ylim([0, ax.YLim(2)])
savefig(gcf, 'figures/follower_duration_init_prob.fig');

delay_1 = delay_e_all_1{1, 14};
delay_2 = delay_e_all_1{11, 14};
index_1 = index_e_all_1{1, 14};
index_2 = index_e_all_1{11, 14};

M = [(1:length(index_1))', delay_1];
sortM = sortrows(M, 2);
sort_ind_1 = sortM(:, 1);
common_idx_1 = zeros(length(sort_ind_1), 1);
for i = 1:length(sort_ind_1)
        if sum(index_1(sort_ind_1(i)) == index_2)==1
                common_idx_1(i) = 1;
        end
end
figure
set(gcf, 'color', [1 1 1])
plot(delay_1(sort_ind_1)*1e3, 1:length(index_1), '-o', 'LineWidth', 1)
hold on
scatter(delay_1(sort_ind_1(common_idx_1==1))*1e3, find(common_idx_1 == 1), 'filled', 'MarkerFaceColor', 'k')
set(gca, 'YDir', 'reverse')
xlabel('Median delay (ms)', 'FontSize', 15);
ylabel('#Neurons', 'FontSize', 15);
ax = gca;
a = get(ax,'Label');
set(ax,'Label', a,'fontsize',15);
set(ax, 'Box', 'off')
savefig(gcf, 'figures/example_sequence_1.fig');

M = [(1:length(index_2))', delay_2];
sortM = sortrows(M, 2);
sort_ind_2 = sortM(:, 1);
common_idx_2 = zeros(length(sort_ind_2), 1);
for i = 1:length(sort_ind_2)
        if sum(index_2(sort_ind_2(i)) == index_1)==1
                common_idx_2(i) = 1;
        end
end
figure
set(gcf, 'color', [1 1 1])
plot(delay_2(sort_ind_2)*1e3, 1:length(index_2), '-o', 'LineWidth', 1)
hold on
scatter(delay_2(sort_ind_2(common_idx_2==1))*1e3, find(common_idx_2 == 1), 'filled', 'MarkerFaceColor', 'k')
set(gca, 'YDir', 'reverse')
xlabel('Median delay (ms)', 'FontSize', 15);
ylabel('#Neurons', 'FontSize', 15);
ax = gca;
a = get(ax,'Label');
set(ax,'Label', a,'fontsize',15);
set(ax, 'Box', 'off')
savefig(gcf, 'figures/example_sequence_2.fig');

%}
