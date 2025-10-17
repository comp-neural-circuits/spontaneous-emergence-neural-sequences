load('followers_corr_uncorr.mat')
mkdir('figures_followers_corr_uncorr');

figure
scatter(cell2mat(followers_corr), cell2mat(followers_uncorr), 'filled')
hold on
ax = gca;
a = get(ax,'Label');
set(gca,'Label', a,'fontsize',14)
set(gcf, 'color', [1 1 1])
plot([0 10], [0 10]*11, 'k', 'LineWidth', 1)
ylim([0 10])
p_num_follower = signrank(cell2mat(followers_uncorr), 11*cell2mat(followers_corr), 'tail', 'left');
xlabel('Number of followers in the block')
ylabel('Number of followers outside the block')
text(6, 8, ['P = ', num2str(round(p_num_follower, ceil(-log10(p_num_follower))+1)*10^ceil(-log10(p_num_follower))), ...
        ' x 10^{-', num2str(ceil(-log10(p_num_follower))), '}'],'fontsize', 16)
savefig(gcf, 'figures_followers_corr_uncorr/num_fol.fig')

figure
subplot(2,1,1)
histogram(cell2mat(all_prob_corr), 'Normalization', 'pdf', 'BinWidth', 0.05, 'FaceColor', ...
        [0 0 0.7], 'LineStyle', 'none')
ylabel('Probability density')
ax = gca;
a = get(ax,'Label');
set(gca,'Label', a,'fontsize',12)
set(gcf, 'color', [1 1 1])
xlim([0 1.05])
legend('In the block')
subplot(2,1,2)
histogram(cell2mat(all_prob_uncorr), 'Normalization', 'pdf',  'BinWidth', 0.05, 'FaceColor', ...
        [0.3 0.3 0.3], 'LineStyle', 'none')
xlabel('Reliability of followers')
ylabel('Probability density')
set(gcf, 'color', [1 1 1]);
ax = gca;
a = get(ax,'Label');
set(gca,'Label', a,'fontsize',12)
set(gcf, 'color', [1 1 1])
xlim([0 1.05])
legend('Out of the block')
p_value_prob = ranksum(cell2mat(all_prob_corr), cell2mat(all_prob_uncorr), 'tail', 'right');
text(0.5, 2.5, ['P = ', num2str(round(p_value_prob, ceil(-log10(p_value_prob))+1)*10^ceil(-log10(p_value_prob))), ...
        ' x 10^{-', num2str(ceil(-log10(p_value_prob))), '}'], 'FontSize', 15)
savefig(gcf, 'figures_followers_corr_uncorr/prob_fol.fig')

figure
subplot(2,1,1)
histogram(cell2mat(all_delay_corr)*1e3, 'Normalization', 'pdf', 'BinWidth', 5, 'FaceColor', ...
        [0 0 0.7], 'LineStyle', 'none')
ylabel('Probability density')
ax = gca;
a = get(ax,'Label');
set(gca,'Label', a,'fontsize',12)
set(gcf, 'color', [1 1 1])
xlim([0 60])
legend('In the block')
subplot(2,1,2)
histogram(cell2mat(all_delay_uncorr)*1e3, 'Normalization', 'pdf',  'BinWidth', 5, 'FaceColor', ...
        [0.3 0.3 0.3], 'LineStyle', 'none')
xlabel('Median delay of followers (ms)')
ylabel('Probability density')
set(gcf, 'color', [1 1 1]);
ax = gca;
a = get(ax,'Label');
set(gca,'Label', a,'fontsize',12)
set(gcf, 'color', [1 1 1])
xlim([0 60])
legend('Out of the block')
p_value_delay = ranksum(cell2mat(all_delay_corr), cell2mat(all_delay_uncorr), 'tail', 'left');
text(40, 0.03, ['P = ', num2str(round(p_value_delay, ceil(-log10(p_value_delay))+1)*10^ceil(-log10(p_value_delay))), ...
        ' x 10^{-', num2str(ceil(-log10(p_value_delay))), '}'], 'FontSize', 15)
savefig(gcf, 'figures_followers_corr_uncorr/delay_fol.fig')

