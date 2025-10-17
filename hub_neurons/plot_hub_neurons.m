load('2/spikes_count.mat');
load('2/final_net.mat');

if ~isfolder('hub_neurons_figure')
        mkdir('hub_neurons_figure')
end

smooth_win = 100;
ckpt_num1 = 6000;     %total number of checkpoints
re = zeros(length(spikes_e{ckpt_num1, 1}), round(ckpt_num1/smooth_win));
ri = zeros(length(spikes_i{ckpt_num1, 1}), round(ckpt_num1/smooth_win));
for ckpt = 1:round(ckpt_num1/smooth_win)
        for frame = 1:smooth_win
                re(:, ckpt) = re(:, ckpt) + spikes_e{(ckpt-1)*smooth_win+frame, 1}/smooth_win;
                ri(:, ckpt) = ri(:, ckpt) + spikes_i{(ckpt-1)*smooth_win+frame, 1}/smooth_win;
        end
end

EE = retrieve_sp(EE_seg);
IE = retrieve_sp(IE_seg);

figure
scatter(re(:, end), log10(max(EE)), 'filled')
xlabel('Firing rate of excitatory neurons (Hz)')
ylabel('Maximal output E-E weight (nS)')
yticks([1,2])
yticklabels({'10', '100'})
yticks([1, log10(50), 2])
yticklabels({'10', '50', '100'})
ax = gca;
a = get(ax,'Label');
set(gca,'Label', a,'fontsize',15)
set(gcf, 'color', [1 1 1])
savefig(gcf, 'hub_neurons_figure/EE_output.fig')
saveas(gcf, 'hub_neurons_figure/EE_output.pdf')

figure
scatter(re(:, end), log10(max(EE, [], 2)), 'filled', 'b')
xlabel('Firing rate of excitatory neurons (Hz)')
ylabel('Maximal incoming E-E weight (nS)')
yticks([1,2])
yticklabels({'10', '100'})
yticks([1, log10(50), 2])
yticklabels({'10', '50', '100'})
ax = gca;
a = get(ax,'Label');
set(gca,'Label', a,'fontsize',15)
set(gcf, 'color', [1 1 1])
savefig(gcf, 'hub_neurons_figure/EE_incoming.fig')
saveas(gcf, 'hub_neurons_figure/EE_incoming.pdf')

figure
scatter(ri(:, end), max(IE), 'filled', 'r')
xlabel('Firing rate of inhibitory neurons (Hz)')
ylabel('Maximal output I-E weight (nS)')
ylim([15 35])
ax = gca;
a = get(ax,'Label');
set(gca,'Label', a,'fontsize',15)
set(gcf, 'color', [1 1 1])
savefig(gcf, 'hub_neurons_figure/IE_output.fig')
saveas(gcf, 'hub_neurons_figure/IE_output.pdf')
