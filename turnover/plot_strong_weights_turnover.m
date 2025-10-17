load('strong_weights_turnover/strong_weights_all_resolution.mat');

resolution = 1000;

ckpt_num1 = 1001;     %total number of checkpoints

trials = 10;

start_point = 7100;   %starting point of new net
early_interval_w = 100;  %time interval to save weights for the old net
t_offset = 100;
early_t_offset = 100;
offset = 8000;    %time offset

all_ckpt = (ckpt_num1 - 4*resolution):resolution:ckpt_num1;
        
if all_ckpt(1) < 1
        plot_ckpt = (early_t_offset+1-start_point):early_interval_w:ckpt_num1;
else
        plot_ckpt = 1:ckpt_num1;
end

avg_EE_ratio_strong = zeros(length(plot_ckpt), length(all_ckpt));
avg_IE_ratio_strong = zeros(length(plot_ckpt), length(all_ckpt));

for i = 1:length(all_ckpt)
        avg_EE_ratio_strong(:, i) = mean(cell2mat(EE_ratio_strong(i, :)), 2);
        avg_IE_ratio_strong(:, i) = mean(cell2mat(IE_ratio_strong(i, :)), 2);
end
        
fig_legend = cell(length(all_ckpt), 1);

plt_avgEE = zeros(length(all_ckpt), 1);

figure()
set(gcf, 'color', [1 1 1]);
hold on
for i = 1:length(all_ckpt)
        adjust_factor = (i-1)/(length(all_ckpt)-1);
        plt_avgEE(i) = plot(plot_ckpt + start_point - early_t_offset + offset - 1, avg_EE_ratio_strong(:, i), 'color', ...
                [0.8*(1 - adjust_factor), 0.8*(1 - adjust_factor), 1 - 0.5*adjust_factor], 'LineWidth', 2);
        fig_legend{i} = ['t = ', num2str((all_ckpt(i)-1) + start_point - early_t_offset + offset), 's'];
end
xlabel('Time (s)', 'FontSize', 15);
ylabel('Average effective weight (nS/mV)', 'FontSize', 15);
title('Strongest E-E connections', 'FontSize', 15);
legend(plt_avgEE, fig_legend, 'Location', 'northwest')
ax = gca;
a = get(ax,'Label');
set(gca,'Label', a,'fontsize',15);
savefig(gcf, ['strong_weights_turnover/ee_strong_avg_', num2str(resolution), 's.fig']);

plt_avgIE = zeros(length(all_ckpt), 1);
figure()
set(gcf, 'color', [1 1 1]);
for i = 1:length(all_ckpt)
        adjust_factor = (i-1)/(length(all_ckpt)-1);
        plt_avgIE(i) = plot(plot_ckpt + start_point - early_t_offset + offset - 1, avg_IE_ratio_strong(:, i), 'color', ...
                [1 - 0.5*adjust_factor, 0.8*(1 - adjust_factor), 0.8*(1 - adjust_factor)], 'LineWidth', 2);
        hold on
end
xlabel('Time (s)', 'FontSize', 15);
ylabel('Average effective weight (nS/mV)', 'FontSize', 15);
title('Strongest I-E connections', 'FontSize', 15);
legend(plt_avgIE, fig_legend, 'Location', 'northwest')
ax = gca;
a = get(ax,'Label');
set(gca,'Label', a,'fontsize',15);
savefig(gcf, ['strong_weights_turnover/ie_strong_avg_', num2str(resolution), 's.fig']);

save('strong_weights_turnover/strong_weights_turnover_traces.mat', 'avg_EE_ratio_strong', 'avg_IE_ratio_strong');



