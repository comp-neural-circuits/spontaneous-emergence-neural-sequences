if ~isfolder('figures_followers_all_blocks')
        mkdir('figures_followers_all_blocks')
end

load('followers_all_blocks.mat')

seg_num = size(all_delay, 2);

prob_avg = zeros(seg_num, 1);
prob_sem = zeros(seg_num, 1);
prob_std = zeros(seg_num, 1);
prob_min = zeros(seg_num, 1);
prob_max = zeros(seg_num, 1);
delay_avg = zeros(seg_num, 1);
delay_sem = zeros(seg_num, 1);
delay_std = zeros(seg_num, 1);
delay_min = zeros(seg_num, 1);
delay_max = zeros(seg_num, 1);
num_avg = zeros(seg_num, 1);
num_sem = zeros(seg_num, 1);
num_std = zeros(seg_num, 1);
num_min = zeros(seg_num, 1);
num_max = zeros(seg_num, 1);

for i = 1:seg_num
        prob_avg(i) = mean(cell2mat(all_prob(:, i)));
        prob_sem(i) = std(cell2mat(all_prob(:, i)))/sqrt(length(cell2mat(all_prob(:, i))));
        prob_std(i) = std(cell2mat(all_prob(:, i)));
        delay_avg(i) = mean(cell2mat(all_delay(:, i))) * 1e3;
        delay_sem(i) = std(cell2mat(all_delay(:, i)))/sqrt(length(cell2mat(all_delay(:, i)))) * 1e3;
        delay_std(i) = std(cell2mat(all_delay(:, i))) * 1e3;
        num_avg(i) = mean(cell2mat(num_followers(:, i)));
        num_sem(i) = std(cell2mat(num_followers(:, i)))/sqrt(length(cell2mat(num_followers(:, i))));
        num_std(i) = std(cell2mat(num_followers(:, i)));
        num_min(i) = min(cell2mat(num_followers(:, i)));
        num_max(i) = max(cell2mat(num_followers(:, i)));
        if isempty(cell2mat(all_prob(:, i)))
                prob_min(i) = NaN;
                prob_max(i) = NaN;
                delay_min(i) = NaN;
                delay_max(i) = NaN;
        else
                prob_min(i) = min(cell2mat(all_prob(:, i)));
                prob_max(i) = max(cell2mat(all_prob(:, i)));
                delay_min(i) = min(cell2mat(all_delay(:, i))) * 1e3;
                delay_max(i) = max(cell2mat(all_delay(:, i))) * 1e3;
        end
end

figure
set(gcf, 'color', [1 1 1])
fill([(1:9), flip(1:9)]', [prob_avg(1:9) - prob_std(1:9); flipud(prob_avg(1:9) + prob_std(1:9))], ...
        [101, 200, 208]/255, 'LineStyle', 'none')
%fill([(1:9), flip(1:9)]', [prob_min(1:9); flipud(prob_max(1:9))], [101, 200, 208]/255, 'LineStyle', 'none')
hold on
fill([(1:9), flip(1:9)]', [prob_avg(1:9) - prob_sem(1:9); flipud(prob_avg(1:9) + prob_sem(1:9))], ...
        [0.7 0.7 1], 'LineStyle', 'none')
plot(prob_avg, 'b', 'LineWidth', 2)
fill([(11:12), flip(11:12)]', [prob_avg(11:12) - prob_std(11:12); flipud(prob_avg(11:12) + prob_std(11:12))], ...
        [101, 200, 208]/255, 'LineStyle', 'none')
%fill([(11:12), flip(11:12)]', [prob_min(11:12); flipud(prob_max(11:12))], [101, 200, 208]/255, 'LineStyle', 'none')
fill([(11:12), flip(11:12)]', [prob_avg(11:12) - prob_sem(11:12); flipud(prob_avg(11:12) + prob_sem(11:12))], ...
        [0.7 0.7 1], 'LineStyle', 'none')
plot(prob_avg, 'b', 'LineWidth', 2)
legend('Std', 's.e.m.', 'Average')
xlabel('No. of blocks starting from the source neuron')
ylabel('Reliability of followers')
ax = gca;
a = get(ax,'Label');
set(gca,'Label', a,'fontsize',12)
set(gca, 'Box', 'off')
xlim([1 12])
ylim([0 ax.YLim(2)])
savefig('figures_followers_all_blocks/prob_avg_all_blocks.fig')

figure
set(gcf, 'color', [1 1 1])
fill([(1:9), flip(1:9)]', [delay_avg(1:9) - delay_std(1:9); flipud(delay_avg(1:9) + delay_std(1:9))], ...
        [101, 200, 208]/255, 'LineStyle', 'none')
%fill([(1:9), flip(1:9)]', [delay_min(1:9); flipud(delay_max(1:9))], [101, 200, 208]/255, 'LineStyle', 'none')
hold on
fill([(1:9), flip(1:9)]', [delay_avg(1:9) - delay_sem(1:9); flipud(delay_avg(1:9) + delay_sem(1:9))], ...
        [0.7 0.7 1], 'LineStyle', 'none')
plot(delay_avg, 'b', 'LineWidth', 2)
fill([(11:12), flip(11:12)]', [delay_avg(11:12) - delay_std(11:12); flipud(delay_avg(11:12) + delay_std(11:12))], ...
        [101, 200, 208]/255, 'LineStyle', 'none')
%fill([(11:12), flip(11:12)]', [delay_min(11:12); flipud(delay_max(11:12))], [101, 200, 208]/255, 'LineStyle', 'none')
fill([(11:12), flip(11:12)]', [delay_avg(11:12) - delay_sem(11:12); flipud(delay_avg(11:12) + delay_sem(11:12))], ...
        [0.7 0.7 1], 'LineStyle', 'none')
plot(delay_avg, 'b', 'LineWidth', 2)
legend('Std', 's.e.m.', 'Average')
xlabel('No. of blocks starting from the source neuron')
ylabel('Median delay of followers (ms)')
ax = gca;
a = get(ax,'Label');
set(gca,'Label', a,'fontsize',12)
set(gca, 'Box', 'off')
xlim([1 12])
ylim([0 ax.YLim(2)])
savefig('figures_followers_all_blocks/delay_avg_all_blocks.fig')

figure
set(gcf, 'color', [1 1 1])
%fill([(1:9), flip(1:9)]', [num_avg(1:9) - num_std(1:9); flipud(num_avg(1:9) + num_std(1:9))], ...
    %    [101, 200, 208]/255, 'LineStyle', 'none')
%fill([(1:9), flip(1:9)]', [num_min(1:9); flipud(num_max(1:9))], [101, 200, 208]/255, 'LineStyle', 'none')
hold on
fill([(1:9), flip(1:9)]', [num_avg(1:9) - num_sem(1:9); flipud(num_avg(1:9) + num_sem(1:9))], ...
        [0.7 0.7 1], 'LineStyle', 'none')
plot(num_avg, 'b', 'LineWidth', 2)
%fill([(11:12), flip(11:12)]', [num_avg(11:12) - num_std(11:12); flipud(num_avg(11:12) + num_std(11:12))], ...
    %    [101, 200, 208]/255, 'LineStyle', 'none')
%fill([(11:12), flip(11:12)]', [num_min(11:12); flipud(num_max(11:12))], [101, 200, 208]/255, 'LineStyle', 'none')
fill([(11:12), flip(11:12)]', [num_avg(11:12) - num_sem(11:12); flipud(num_avg(11:12) + num_sem(11:12))], ...
        [0.7 0.7 1], 'LineStyle', 'none')
plot(num_avg, 'b', 'LineWidth', 2)
legend('s.e.m.', 'Average')
xlabel('No. of blocks starting from the source neuron')
ylabel('Number of followers')
ax = gca;
a = get(ax,'Label');
set(gca,'Label', a,'fontsize',12)
set(gca, 'Box', 'off')
xlim([1 12])
ylim([0 ax.YLim(2)])
savefig('figures_followers_all_blocks/num_avg_all_blocks.fig')




