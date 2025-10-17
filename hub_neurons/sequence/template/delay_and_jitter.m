%first load the sequence data

if ~isfolder('sequence')
        mkdir('sequence');
end

Ne = 1200;
Ni = 240;
sim_num = 10;
dt = 1e-4;

T = 43;
frames = round(ceil(T)/dt);

index_e = 1:Ne;
index_i = 1:Ni;
for i = 1:sim_num
        load(['sequence_', num2str(i), '.mat'])
        index_e = intersect(index_e, follower_index_e);
        index_i = intersect(index_i, follower_index_i);
end

stim_num = 100;

delay_e = cell(sim_num, 1);
delay_i = cell(sim_num, 1);

bin_width = 0.1;

parfor i = 1:sim_num
        
        disp(i);
        
        delay_e{i} = zeros(length(index_e), stim_num);
        delay_i{i} = zeros(length(index_i), stim_num);
        
        folder = num2str(i);

        xe_seq = zeros(length(index_e), frames);
        xi_seq = zeros(length(index_i), frames);

        stim_time = stim_onset + (0: (stim_num - 1))*stim_period;

        for seg_index = 1:ceil(T/bin_width)
                frame = (1+(seg_index-1)*round(bin_width/dt)):(seg_index*round(bin_width/dt));
                xe_seg = h5read_sparsematrix([folder, '/all_x.h5'], [num2str(seg_index), '/xe']);
                xe = retrieve_sp(xe_seg);
                xe_seq(:, frame) = xe(index_e, 1:end-1);
                xi_seg = h5read_sparsematrix([folder, '/all_x.h5'], [num2str(seg_index), '/xi']);
                xi = retrieve_sp(xi_seg);
                xi_seq(:, frame) = xi(index_i, 1:end-1);
        end

        for neuron = 1:length(index_e)
                spike_time = find(xe_seq(neuron, :) == 1) * dt;
                for trial = 1:stim_num
                        follower_time = min(spike_time(spike_time > stim_time(trial)));
                        if follower_time < stim_time(trial) + effect
                                delay_e{i}(neuron, trial) = follower_time - stim_time(trial);
                        else
                                delay_e{i}(neuron, trial) = -1;
                        end
                end
        end
        
         for neuron = 1:length(index_i)
                spike_time = find(xi_seq(neuron, :) == 1) * dt;
                for trial = 1:stim_num
                        follower_time = min(spike_time(spike_time > stim_time(trial)));
                        if follower_time < stim_time(trial) + effect
                                delay_i{i}(neuron, trial) = follower_time - stim_time(trial);
                        else
                                delay_i{i}(neuron, trial) = -1;
                        end
                end
        end
end
                        
delay_e_all = cell2mat(delay_e');
delay_i_all = cell2mat(delay_i');

med_delay_e = zeros(length(index_e), 1);
med_delay_i = zeros(length(index_i), 1);
jitter_e = zeros(length(index_e), 1);
jitter_i = zeros(length(index_i), 1);

for neuron = 1:length(index_e)
        med_delay_e(neuron) = median(delay_e_all(neuron, delay_e_all(neuron, :) > 0));
        jitter_e(neuron) = std(delay_e_all(neuron, delay_e_all(neuron, :) > 0));
end
for neuron = 1:length(index_i)
        med_delay_i(neuron) = median(delay_i_all(neuron, delay_i_all(neuron, :) > 0));
        jitter_i(neuron) = std(delay_i_all(neuron, delay_i_all(neuron, :) > 0));
end

prob_follower_e = mean(delay_e_all > 0, 2);
prob_follower_i = mean(delay_i_all > 0, 2);

sort_delay_e = [med_delay_e, (1:length(index_e))'];
sort_delay_e = sortrows(sort_delay_e, 1);
sorted_idx_e = sort_delay_e(:, 2);

figure;
fill([med_delay_e(sorted_idx_e) - jitter_e(sorted_idx_e); ...
        flip(med_delay_e(sorted_idx_e) + jitter_e(sorted_idx_e))]*1e3, ...
        [(1:length(index_e))'; flip(1:length(index_e))'], 'c', 'LineStyle', 'none');
hold on
plot(med_delay_e(sorted_idx_e)*1e3, 1:length(index_e), 'b', 'LineWidth', 2)
set(gca, 'YDir', 'reverse')
xlabel('Median delay (ms)', 'FontSize', 15);
ylabel('#Neurons', 'FontSize', 15);
set(gcf, 'color', [1 1 1])
ax = gca;
a = get(ax,'Label');
set(gca,'Label', a,'fontsize',15);
title('Excitatory', 'FontSize', 15)
if length(index_e) > 1
        ylim([1 length(index_e)])
end
savefig(gcf, 'sequence/e_sequence.fig');

sort_delay_i = [med_delay_i, (1:length(index_i))'];
sort_delay_i = sortrows(sort_delay_i, 1);
sorted_idx_i = sort_delay_i(:, 2);

fig = figure;
fill([med_delay_i(sorted_idx_i) - jitter_i(sorted_idx_i); ...
        flip(med_delay_i(sorted_idx_i) + jitter_i(sorted_idx_i))]*1e3, ...
        [(1:length(index_i))'; flip(1:length(index_i))'], 'm', 'LineStyle', 'none');
hold on
plot(med_delay_i(sorted_idx_i)*1e3, 1:length(index_i), 'r', 'LineWidth', 2)
set(gca, 'YDir', 'reverse')
xlabel('Median delay (ms)', 'FontSize', 15);
ylabel('#Neurons', 'FontSize', 15);
set(gcf, 'color', [1 1 1])
ax = gca;
a = get(ax,'Label');
set(gca,'Label', a,'fontsize',15);
title('Inhibitory', 'FontSize', 15)
if length(index_i) > 1
        ylim([1 length(index_i)])
end
savefig(gcf, 'sequence/i_sequence.fig');

figure
set(gcf, 'color', [1 1 1])
scatter(prob_follower_e, med_delay_e*1e3, 'filled', 'MarkerFaceColor', 'b')
hold on
scatter(prob_follower_i, med_delay_i*1e3, 'filled', 'MarkerFaceColor', 'r')
xlabel('Probability of being a follower', 'FontSize', 15);
ylabel('Median delay (ms)', 'FontSize', 15);
h = legend('Excitatory', 'Inhibitory');
set(h, 'FontSize', 15)
ax = gca;
a = get(ax,'Label');
set(gca,'Label', a,'fontsize',15);
savefig('sequence/prob_delay.fig');

figure
set(gcf, 'color', [1 1 1])
scatter(prob_follower_e, jitter_e*1e3, 'filled', 'MarkerFaceColor', 'b')
hold on
scatter(prob_follower_i, jitter_i*1e3, 'filled', 'MarkerFaceColor', 'r')
xlabel('Probability of being a follower', 'FontSize', 15);
ylabel('Jitter (ms)', 'FontSize', 15);
h = legend('Excitatory', 'Inhibitory');
set(h, 'FontSize', 15)
ax = gca;
a = get(ax,'Label');
set(gca,'Label', a,'fontsize',15);
savefig('sequence/prob_jitter.fig');

figure
set(gcf, 'color', [1 1 1])
scatter(med_delay_e*1e3, jitter_e*1e3, 'filled', 'MarkerFaceColor', 'b')
hold on
scatter(med_delay_i*1e3, jitter_i*1e3, 'filled', 'MarkerFaceColor', 'r')
xlabel('Median delay (ms)', 'FontSize', 15);
ylabel('Jitter (ms)', 'FontSize', 15);
h = legend('Excitatory', 'Inhibitory', 'Location', 'southeast');
set(h, 'FontSize', 15)
ax = gca;
a = get(ax,'Label');
set(gca,'Label', a,'fontsize',15);
savefig('sequence/delay_jitter.fig');
        
save('sequence/delay_and_jitter.mat', 'index_e', 'index_i', 'sorted_idx_e', 'sorted_idx_i', ...
        'med_delay_e', 'med_delay_i', 'jitter_e', 'jitter_i', 'prob_follower_e', 'prob_follower_i');


        
        
