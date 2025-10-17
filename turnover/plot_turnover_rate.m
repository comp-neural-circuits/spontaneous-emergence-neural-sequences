load('strong_weights_turnover/turnover_rate.mat');

%cut data
index = 1:12;
resolution = resolution(index);
n0_ee = n0_ee(:, index);
n1_ee = n1_ee(:, index);
n0_ie = n0_ie(:, index);
n1_ie = n1_ie(:, index);
ratio_n0_n1 = ratio_n0_n1(index, :, :, :);

trials = size(ratio_n0_n1, 2);

%turnover rate

remain_rate_ee = zeros(length(resolution), trials);
remain_rate_ie = zeros(length(resolution), trials);

figure
set(gcf, 'color', [1 1 1])
for i = 1:length(resolution)
        remain_rate_ee(i, :) =  sum(cell2mat(n0_ee(:, i)'), 1) ./ sum(cell2mat(n0_ee(:, i)') + cell2mat(n1_ee(:, i)'), 1);
        remain_rate_ie(i, :) =  sum(cell2mat(n0_ie(:, i)'), 1) ./ sum(cell2mat(n0_ie(:, i)') + cell2mat(n1_ie(:, i)'), 1);
        scatter(resolution(i)*ones(trials, 1), remain_rate_ee(i, :), 'b', 'filled')
        hold on
        scatter(resolution(i)*ones(trials, 1), remain_rate_ie(i, :), 'r', 'filled')
end

fit_decay_ee = fit(repmat(resolution, [1, trials])', remain_rate_ee(:), 'exp(-x/a)', 'Start', 1000);
fit_decay_ie = fit(repmat(resolution, [1, trials])', remain_rate_ie(:), 'exp(-x/a)',  'Start', 1000);
plot_x = 1:max(resolution);
plot(plot_x, exp(-plot_x/fit_decay_ee.a), 'k', 'LineWidth', 1);
plot(plot_x, exp(-plot_x/fit_decay_ie.a), 'k', 'LineWidth', 1);
legend('E-E', 'I-E')
xlabel('Interval (s)')
ylabel('Ratio of remaining strong connections')
ax = gca;
a = get(ax,'Label');
set(gca,'Label', a,'fontsize',15)
text(1200, 0.6, ['\tau_{ee} = ', num2str(round(fit_decay_ee.a)), 's'], 'FontSize', 15, 'color', 'b')
text(200, 0.3, ['\tau_{ie} = ', num2str(round(fit_decay_ie.a)), 's'], 'FontSize', 15, 'color', 'r')
ylim([0 1])

savefig(gcf, 'strong_weights_turnover/turnover_rate.fig')

