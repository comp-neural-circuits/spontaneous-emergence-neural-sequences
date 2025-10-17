result_dir = 'triplet_turnover/';

if ~isfolder(result_dir)
        mkdir(result_dir);
end

load('triplet_motif_turnover.mat');

trials = size(n0, 1);

motif_types = 5;

%turnover rate

remain_rate = zeros(length(resolution), trials, motif_types);

for motif = 1:motif_types
        for i = 1:length(resolution)
                for trial = 1:trials
                        remain_rate(i, trial, motif) =  sum(n0{trial, i}(:, motif))./(sum(n0{trial, i}(:, motif)) + sum(n1{trial, i}(:, motif)));
                end
        end
        figure
        set(gcf, 'color', [1 1 1])
        scatter(repmat(resolution, [1, trials])', reshape(remain_rate(:, :, motif), trials*length(resolution), 1), 'b', 'filled')
        if motif ~= 5
                hold on
                fit_decay = fit(repmat(resolution, [1, trials])', reshape(remain_rate(:, :, motif), trials*length(resolution), 1)...
                        , 'exp(-x/a)', 'Start', 1000);
                plot_x = 1:max(resolution);
                plot(plot_x, exp(-plot_x/fit_decay.a), 'k', 'LineWidth', 1);
                text(600, 0.6, ['\tau_{motif} = ', num2str(round(fit_decay.a)), 's'], 'FontSize', 15, 'color', 'b')
        end
         xlabel('Time (s)')
         ylabel('Ratio of remaining motifs')
         ax = gca;
         a = get(ax,'Label');
         set(gca,'Label', a,'fontsize',15)
         savefig(gcf, [result_dir, 'turnover_rate_tri_motif_', num2str(motif), '.fig'])
end

save('decay_ratio.mat', 'resolution', 'remain_rate')
