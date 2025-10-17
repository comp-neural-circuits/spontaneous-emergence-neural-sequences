%Fitting: exp(-t/tau)

load('strong_weights_turnover/turnover_rate.mat')
figure
hold on
set(gcf, 'color', [1 1 1])

%cut data
index = 1:6;
n0_ee = n0_ee(:, index);
n1_ee = n1_ee(:, index);
resolution = resolution(index);

trials = size(n0_ee, 1);
remain_rate_ee = zeros(length(resolution), trials);

for i = 1:length(resolution)
        remain_rate_ee(i, :) =  sum(cell2mat(n0_ee(:, i)'), 1)./sum(cell2mat(n0_ee(:, i)') + cell2mat(n1_ee(:, i)'), 1);
        scatter(resolution(i)*ones(trials, 1), remain_rate_ee(i, :), 'filled', 'MarkerFaceColor', [117, 76, 36]/255)
end

fit_decay_ee = fit(repmat(resolution, [1, trials])', remain_rate_ee(:), 'exp(-x/a)', 'Start', 1000);
plot_x = 1:max(resolution);
plot(plot_x, exp(-plot_x/fit_decay_ee.a), 'color', [117, 76, 36]/255, 'LineWidth', 1);
xlabel('Interval (s)')
ylabel('Ratio of remaining strong connections')
ax = gca;
a = get(ax,'Label');
set(gca,'Label', a,'fontsize',15)
text(200, 0.6, ['\tau_{ee} = ', num2str(round(fit_decay_ee.a)), 's'], 'FontSize', 15, 'color', [117, 76, 36]/255)
ylim([0 1])

%%
load('strong_weights_turnover/turnover_corr_rate.mat')
set(gcf, 'color', [1 1 1])

%cut data
index = 1:6;
n0_ee = n0_ee(:, index);
n1_ee = n1_ee(:, index);
resolution = resolution(index);

trials = size(n0_ee, 1);
remain_rate_ee = zeros(length(resolution), trials);

for i = 1:length(resolution)
        remain_rate_ee(i, :) =  sum(cell2mat(n0_ee(:, i)'), 1)./sum(cell2mat(n0_ee(:, i)') + cell2mat(n1_ee(:, i)'), 1);
        scatter(resolution(i)*ones(trials, 1), remain_rate_ee(i, :), 'filled', 'MarkerFaceColor', [46, 49, 146]/255)
end

fit_decay_ee = fit(repmat(resolution, [1, trials])', remain_rate_ee(:), 'exp(-x/a)', 'Start', 1000);
plot_x = 1:max(resolution);
plot(plot_x, exp(-plot_x/fit_decay_ee.a), 'color', [46, 49, 146]/255, 'LineWidth', 1);
xlabel('Interval (s)')
ylabel('Ratio of remaining strong connections')
ax = gca;
a = get(ax,'Label');
set(gca,'Label', a,'fontsize',15)
text(200, 0.5, ['\tau_{ee} = ', num2str(round(fit_decay_ee.a)), 's'], 'FontSize', 15, 'color', [46, 49, 146]/255)
ylim([0 1])

%%
load('../turnover/strong_weights_turnover/turnover_rate.mat')
set(gcf, 'color', [1 1 1])

%cut data
index = 7:12;
n0_ee = n0_ee(:, index);
n1_ee = n1_ee(:, index);
resolution = resolution(index);

trials = size(n0_ee, 1);
remain_rate_ee = zeros(length(resolution), trials);

for i = 1:length(resolution)
        remain_rate_ee(i, :) =  sum(cell2mat(n0_ee(:, i)'), 1)./sum(cell2mat(n0_ee(:, i)') + cell2mat(n1_ee(:, i)'), 1);
        scatter(resolution(i)*ones(trials, 1), remain_rate_ee(i, :), 'filled', 'MarkerFaceColor', 'k')
end

fit_decay_ee = fit(repmat(resolution, [1, trials])', remain_rate_ee(:), 'exp(-x/a)', 'Start', 1000);
plot_x = 1:max(resolution);
plot(plot_x, exp(-plot_x/fit_decay_ee.a), 'color', 'k', 'LineWidth', 1);
xlabel('Interval (s)')
ylabel('Ratio of remaining strong connections')
ax = gca;
a = get(ax,'Label');
set(gca,'Label', a,'fontsize',15)
text(200, 0.4, ['\tau_{ee} = ', num2str(round(fit_decay_ee.a)), 's'], 'FontSize', 15, 'color', 'k')
ylim([0 1])


%%
load('strong_weights_turnover/turnover_uncorr_rate.mat')
set(gcf, 'color', [1 1 1])

%cut data
index = 1:6;
n0_ee = n0_ee(:, index);
n1_ee = n1_ee(:, index);
resolution = resolution(index);

trials = size(n0_ee, 1);
remain_rate_ee = zeros(length(resolution), trials);

for i = 1:length(resolution)
        remain_rate_ee(i, :) =  sum(cell2mat(n0_ee(:, i)'), 1)./sum(cell2mat(n0_ee(:, i)') + cell2mat(n1_ee(:, i)'), 1);
        scatter(resolution(i)*ones(trials, 1), remain_rate_ee(i, :), 'filled', 'MarkerFaceColor', [1 1 1]*0.5)
end

fit_decay_ee = fit(repmat(resolution, [1, trials])', remain_rate_ee(:), 'exp(-x/a)', 'Start', 1000);
plot_x = 1:max(resolution);
plot(plot_x, exp(-plot_x/fit_decay_ee.a), 'color', [0.5, 0.5, 0.5], 'LineWidth', 1);
xlabel('Interval (s)')
ylabel('Ratio of remaining strong connections')
ax = gca;
a = get(ax,'Label');
set(gca,'Label', a,'fontsize',15)
text(200, 0.3, ['\tau_{ee} = ', num2str(round(fit_decay_ee.a)), 's'], 'FontSize', 15, 'color', [0.5, 0.5, 0.5])
ylim([0 1])

%%
load('strong_weights_turnover/turnover_corr_block_rate.mat')
set(gcf, 'color', [1 1 1])

%cut data
index = 1:6;
n0_ee = n0_ee(:, index);
n1_ee = n1_ee(:, index);
resolution = resolution(index);

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

fit_decay_ee = fit(repmat(resolution, [1, trials])', remain_rate_ee(:), 'exp(-x/a)', 'Start', 1000);
plot_x = 1:max(resolution);
plot(plot_x, exp(-plot_x/fit_decay_ee.a), 'color', [237, 30, 121]/255, 'LineWidth', 1);
xlabel('Interval (s)')
ylabel('Ratio of remaining strong connections')
ax = gca;
a = get(ax,'Label');
set(gca,'Label', a,'fontsize',15)
text(200, 0.2, ['\tau_{ee} = ', num2str(round(fit_decay_ee.a)), 's'], 'FontSize', 15, 'color', [237, 30, 121]/255)
ylim([0 1])

savefig(gcf, 'strong_weights_turnover/turnover_all_corr_uncorr.fig')





