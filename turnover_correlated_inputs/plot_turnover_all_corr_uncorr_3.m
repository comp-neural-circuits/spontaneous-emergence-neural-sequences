%Fitting: (x/a+1)^(-b)

load('strong_weights_turnover/turnover_rate.mat')
figure
hold on
set(gcf, 'color', [1 1 1])

trials = size(n0_ee, 1);
remain_rate_ee = zeros(length(resolution), trials);

for i = 1:length(resolution)
        remain_rate_ee(i, :) =  sum(cell2mat(n0_ee(:, i)'), 1)./sum(cell2mat(n0_ee(:, i)') + cell2mat(n1_ee(:, i)'), 1);
        scatter(resolution(i)*ones(trials, 1), remain_rate_ee(i, :), 'filled', 'MarkerFaceColor', [117, 76, 36]/255)
end

fit_decay_ee = fit(repmat(resolution, [1, trials])', remain_rate_ee(:), '(x/a+1)^(-b)', 'Start', [1000, 1], 'Lower', [0, 0]);
plot_x = 1:max(resolution);
plot(plot_x, (plot_x/fit_decay_ee.a+1).^(-fit_decay_ee.b), 'color', [117, 76, 36]/255, 'LineWidth', 1);
xlabel('Interval (s)')
ylabel('Ratio of remaining strong connections')
ax = gca;
a = get(ax,'Label');
set(gca,'Label', a,'fontsize',15)
text(200, 0.6, ['\delta = ', num2str(round(fit_decay_ee.b, 3))], 'FontSize', 15, 'color', [117, 76, 36]/255)
ylim([0 1])

%%
load('strong_weights_turnover/turnover_corr_rate.mat')
set(gcf, 'color', [1 1 1])

trials = size(n0_ee, 1);
remain_rate_ee = zeros(length(resolution), trials);

for i = 1:length(resolution)
        remain_rate_ee(i, :) =  sum(cell2mat(n0_ee(:, i)'), 1)./sum(cell2mat(n0_ee(:, i)') + cell2mat(n1_ee(:, i)'), 1);
        scatter(resolution(i)*ones(trials, 1), remain_rate_ee(i, :), 'filled', 'MarkerFaceColor', [46, 49, 146]/255)
end

fit_decay_ee = fit(repmat(resolution, [1, trials])', remain_rate_ee(:), '(x/a+1)^(-b)', 'Start', [1000, 1], 'Lower', [0, 0]);
plot_x = 1:max(resolution);
plot(plot_x, (plot_x/fit_decay_ee.a+1).^(-fit_decay_ee.b), 'color', [46, 49, 146]/255, 'LineWidth', 1);
xlabel('Interval (s)')
ylabel('Ratio of remaining strong connections')
ax = gca;
a = get(ax,'Label');
set(gca,'Label', a,'fontsize',15)
text(200, 0.5, ['\delta = ', num2str(round(fit_decay_ee.b, 3))], 'FontSize', 15, 'color', [46, 49, 146]/255)
ylim([0 1])

%%
load('../turnover/strong_weights_turnover/turnover_rate.mat')
set(gcf, 'color', [1 1 1])

trials = size(n0_ee, 1);
remain_rate_ee = zeros(length(resolution), trials);

for i = 1:length(resolution)
        remain_rate_ee(i, :) =  sum(cell2mat(n0_ee(:, i)'), 1)./sum(cell2mat(n0_ee(:, i)') + cell2mat(n1_ee(:, i)'), 1);
        scatter(resolution(i)*ones(trials, 1), remain_rate_ee(i, :), 'filled', 'MarkerFaceColor', 'k')
end

fit_decay_ee = fit(repmat(resolution, [1, trials])', remain_rate_ee(:), '(x/a+1)^(-b)', 'Start', [1000, 1], 'Lower', [0, 0]);
plot_x = 1:max(resolution);
plot(plot_x, (plot_x/fit_decay_ee.a+1).^(-fit_decay_ee.b), 'color', 'k', 'LineWidth', 1);
xlabel('Interval (s)')
ylabel('Ratio of remaining strong connections')
ax = gca;
a = get(ax,'Label');
set(gca,'Label', a,'fontsize',15)
text(200, 0.4, ['\delta = ', num2str(round(fit_decay_ee.b, 3))], 'FontSize', 15, 'color', 'k')
ylim([0 1])


%%
load('strong_weights_turnover/turnover_uncorr_rate.mat')
set(gcf, 'color', [1 1 1])

trials = size(n0_ee, 1);
remain_rate_ee = zeros(length(resolution), trials);

for i = 1:length(resolution)
        remain_rate_ee(i, :) =  sum(cell2mat(n0_ee(:, i)'), 1)./sum(cell2mat(n0_ee(:, i)') + cell2mat(n1_ee(:, i)'), 1);
        scatter(resolution(i)*ones(trials, 1), remain_rate_ee(i, :), 'filled', 'MarkerFaceColor', [1 1 1]*0.5)
end

fit_decay_ee = fit(repmat(resolution, [1, trials])', remain_rate_ee(:), '(x/a+1)^(-b)', 'Start', [1000, 1], 'Lower', [0, 0]);
plot_x = 1:max(resolution);
plot(plot_x, (plot_x/fit_decay_ee.a+1).^(-fit_decay_ee.b), 'color', [0.5, 0.5, 0.5], 'LineWidth', 1);
xlabel('Interval (s)')
ylabel('Ratio of remaining strong connections')
ax = gca;
a = get(ax,'Label');
set(gca,'Label', a,'fontsize',15)
text(200, 0.3, ['\delta = ', num2str(round(fit_decay_ee.b, 3))], 'FontSize', 15, 'color', [0.5, 0.5, 0.5])
ylim([0 1])

%%
load('strong_weights_turnover/turnover_corr_block_rate.mat')
set(gcf, 'color', [1 1 1])

n0_ee_avg = cell(size(n0_ee));
n1_ee_avg = cell(size(n1_ee));
for i = 1:length(resolution)
        for trial = 1:trials
                n0_ee_avg{trial, i} = mean(n0_ee{trial, i}, 2);
                n1_ee_avg{trial, i} = mean(n1_ee{trial, i}, 2);
        end
end

trials = size(n0_ee_avg, 1);
remain_rate_ee = zeros(length(resolution), trials);

for i = 1:length(resolution)
        remain_rate_ee(i, :) =  sum(cell2mat(n0_ee_avg(:, i)'), 1) ./ sum(cell2mat(n0_ee_avg(:, i)') + cell2mat(n1_ee_avg(:, i)'), 1);
        scatter(resolution(i)*ones(trials, 1), remain_rate_ee(i, :), 'filled', 'MarkerFaceColor', [237, 30, 121]/255)
end

fit_decay_ee = fit(repmat(resolution, [1, trials])', remain_rate_ee(:), '(x/a+1)^(-b)', 'Start', [1000, 1], 'Lower', [0, 0]);
plot_x = 1:max(resolution);
plot(plot_x, (plot_x/fit_decay_ee.a+1).^(-fit_decay_ee.b), 'color', [237, 30, 121]/255, 'LineWidth', 1);
xlabel('Interval (s)')
ylabel('Ratio of remaining strong connections')
ax = gca;
a = get(ax,'Label');
set(gca,'Label', a,'fontsize',15)
text(200, 0.2, ['\delta = ', num2str(round(fit_decay_ee.b, 3))], 'FontSize', 15, 'color', [237, 30, 121]/255)
ylim([0 1])

savefig(gcf, 'strong_weights_turnover/turnover_all_corr_uncorr_3.fig')





