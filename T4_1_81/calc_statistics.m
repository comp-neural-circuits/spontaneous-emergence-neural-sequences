load('followers.mat');

EE_dir = '../../shifting_inputs/2/1/';
ckpt = 81;

V0 = -0.063;
EE_seg = h5read_sparsematrix([EE_dir, '/all_w.h5'], [num2str(ckpt), '/EE']);
EE = retrieve_sp(EE_seg);
ET = h5read([EE_dir, '/all_w.h5'], ['/', num2str(ckpt), '/ET']);
effw = max(EE./((ET - V0)*ones(1,size(EE, 2)))/1e3);

groups = size(followers_e, 1);
group_width = 0.5;

prob_e = zeros(groups, 1);
prob_i = zeros(groups, 1);
effw_prob = zeros(groups, 1);

for group = 1:groups
        effw_start = (group-1)*group_width;
        effw_end = group*group_width;
        if group == groups
                effw_prob(group) = mean(effw >= effw_start);
        else
                effw_prob(group) = mean((effw >= effw_start) & (effw < effw_end));
        end
        prob_e(group) = mean(followers_e{group, 1} > 0);
        prob_i(group) = mean(followers_i{group, 1} > 0);
end

overall_prob_e = prob_e(effw_prob > 0)' * effw_prob(effw_prob > 0);
overall_prob_i = prob_i(effw_prob > 0)' * effw_prob(effw_prob > 0);

figure;
set(gcf, 'color', [1 1 1])
histogram(effw, 'BinWidth', 0.5, 'FaceColor', 'k')
hold on
yyaxis right
set(gca, 'ycolor', 'b')
plot(linspace(group_width/2, group_width*(groups-1/2), groups), prob_e, 'b', 'LineWidth', 2)
hold on
plot(linspace(group_width/2, group_width*(groups-1/2), groups), prob_i, 'r', 'LineWidth', 2, 'LineStyle', '-')
xlabel('Maximal weight-thresh ratio (nS/mV)')
ylabel('Probability to generate a sequence')
yyaxis left
ylabel('Count of neurons')
ax = gca;
a = get(ax,'Label');
set(gca,'Label', a,'fontsize',15)
        