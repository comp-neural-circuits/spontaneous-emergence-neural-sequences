load('followers_corr_uncorr2.mat')
if ~isfolder('figures_followers_corr_uncorr2')
        mkdir('figures_followers_corr_uncorr2');
end

figure
scatter(cell2mat(followers_corr_diag), cell2mat(followers_uncorr), 'filled')
hold on
scatter(cell2mat(followers_corr_subdiag), cell2mat(followers_uncorr), 'filled')
ax = gca;
a = get(ax,'Label');
set(gca,'Label', a,'fontsize',14)
set(gcf, 'color', [1 1 1])
plot([0 10], [0 10]*10, 'k', 'LineWidth', 1)
ylim([0 10])
p_num_follower = signrank(cell2mat(followers_uncorr), 10*cell2mat(followers_corr_diag), 'tail', 'left');
xlabel('Number of followers in the block')
ylabel('Number of followers outside the block')
text(6, 8, ['P = ', num2str(round(p_num_follower, ceil(-log10(p_num_follower))+1)*10^ceil(-log10(p_num_follower))), ...
        ' x 10^{-', num2str(ceil(-log10(p_num_follower))), '}'],'fontsize', 16)
savefig(gcf, 'figures_followers_corr_uncorr2/num_fol.fig')

figure
subplot(3,1,1)
histogram(cell2mat(all_prob_corr_diag), 'Normalization', 'pdf', 'BinWidth', 0.05, 'FaceColor', ...
        [0 0 0.7], 'LineStyle', 'none')
%ylabel('Probability density')
ax = gca;
a = get(ax,'Label');
set(gca,'Label', a,'fontsize',12)
set(gcf, 'color', [1 1 1])
xlim([0 1.05])
legend('the same block')
subplot(3,1,2)
histogram(cell2mat(all_prob_corr_subdiag), 'Normalization', 'pdf', 'BinWidth', 0.05, 'FaceColor', ...
        [0.7 0 0], 'LineStyle', 'none')
ylabel('Probability density')
ax = gca;
a = get(ax,'Label');
set(gca,'Label', a,'fontsize',12)
set(gcf, 'color', [1 1 1])
xlim([0 1.05])
legend('the next block')
subplot(3,1,3)
histogram(cell2mat(all_prob_uncorr), 'Normalization', 'pdf',  'BinWidth', 0.05, 'FaceColor', ...
        [0.3 0.3 0.3], 'LineStyle', 'none')
xlabel('Reliability of followers')
%ylabel('Probability density')
set(gcf, 'color', [1 1 1]);
ax = gca;
a = get(ax,'Label');
set(gca,'Label', a,'fontsize',12)
set(gcf, 'color', [1 1 1])
xlim([0 1.05])
legend('other blocks')
p_value_prob_diag = ranksum(cell2mat(all_prob_corr_diag), cell2mat(all_prob_uncorr), 'tail', 'right');
subplot(3, 1, 1)
text(0.2, 5, ['P = ', num2str(round(p_value_prob_diag, ceil(-log10(p_value_prob_diag))+1)*10^ceil(-log10(p_value_prob_diag))), ...
        ' x 10^{-', num2str(ceil(-log10(p_value_prob_diag))), '}'], 'FontSize', 15)
p_value_prob_subdiag = ranksum(cell2mat(all_prob_corr_subdiag), cell2mat(all_prob_uncorr), 'tail', 'right');
subplot(3, 1, 2)
text(0.2, 5, ['P = ', num2str(round(p_value_prob_subdiag, ceil(-log10(p_value_prob_subdiag))+1)*10^ceil(-log10(p_value_prob_subdiag))), ...
        ' x 10^{-', num2str(ceil(-log10(p_value_prob_subdiag))), '}'], 'FontSize', 15)
savefig(gcf, 'figures_followers_corr_uncorr2/prob_fol.fig')

figure
subplot(3,1,1)
histogram(cell2mat(all_delay_corr_diag)*1e3, 'Normalization', 'pdf', 'BinWidth', 5, 'FaceColor', ...
        [0 0 0.7], 'LineStyle', 'none')
%ylabel('delayability density')
ax = gca;
a = get(ax,'Label');
set(gca,'Label', a,'fontsize',12)
set(gcf, 'color', [1 1 1])
xlim([0 50])
legend('the same block')
subplot(3,1,2)
histogram(cell2mat(all_delay_corr_subdiag)*1e3, 'Normalization', 'pdf', 'BinWidth', 5, 'FaceColor', ...
        [0.7 0 0], 'LineStyle', 'none')
ylabel('Probability density')
ax = gca;
a = get(ax,'Label');
set(gca,'Label', a,'fontsize',12)
set(gcf, 'color', [1 1 1])
xlim([0 50])
legend('the next block')
subplot(3,1,3)
histogram(cell2mat(all_delay_uncorr)*1e3, 'Normalization', 'pdf',  'BinWidth', 5, 'FaceColor', ...
        [0.3 0.3 0.3], 'LineStyle', 'none')
xlabel('Median delay of followers')
%ylabel('delayability density')
set(gcf, 'color', [1 1 1]);
ax = gca;
a = get(ax,'Label');
set(gca,'Label', a,'fontsize',12)
set(gcf, 'color', [1 1 1])
xlim([0 70])
legend('other blocks')
p_value_delay_diag = ranksum(cell2mat(all_delay_corr_diag), cell2mat(all_delay_uncorr), 'tail', 'left');
subplot(3, 1, 1)
text(25, 0.03, ['P = ', num2str(round(p_value_delay_diag, ceil(-log10(p_value_delay_diag))+1)*10^ceil(-log10(p_value_delay_diag))), ...
        ' x 10^{-', num2str(ceil(-log10(p_value_delay_diag))), '}'], 'FontSize', 15)
p_value_delay_subdiag = ranksum(cell2mat(all_delay_corr_subdiag), cell2mat(all_delay_uncorr), 'tail', 'left');
subplot(3, 1, 2)
text(30, 0.025, ['P = ', num2str(round(p_value_delay_subdiag, ceil(-log10(p_value_delay_subdiag))+1)*10^ceil(-log10(p_value_delay_subdiag))), ...
        ' x 10^{-', num2str(ceil(-log10(p_value_delay_subdiag))), '}'], 'FontSize', 15)
savefig(gcf, 'figures_followers_corr_uncorr2/delay_fol.fig')

