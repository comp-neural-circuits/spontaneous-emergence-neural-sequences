%%
%visualization

load('sequence/delay_and_jitter.mat')

dt = 1e-4;
folder = '2';

samples = 10;

xe_all = zeros(length(index_e), round(0.1/dt));
xi_all = zeros(length(index_i), round(0.1/dt));
stim_onset = 3.1;
bin_width = 0.1;
for i = 1:round(0.1/bin_width)
        file_num = i+round(stim_onset/bin_width);
        xe_seg = h5read_sparsematrix([folder, '/all_x.h5'], [num2str(file_num), '/xe']);
        xe = retrieve_sp(xe_seg);
        xi_seg = h5read_sparsematrix([folder, '/all_x.h5'], [num2str(file_num), '/xi']);
        xi = retrieve_sp(xi_seg);
        xe_all(:, (1+round((i-1)*bin_width/dt)):round(i*bin_width/dt)) = xe(index_e, 1:end-1);
        xi_all(:, (1+round((i-1)*bin_width/dt)):round(i*bin_width/dt)) = xi(index_i, 1:end-1);
end

first_spike = zeros(length(index_e), 1);
for i = 1:length(index_e)
        if sum(xe(index_e(i), :)) > 0
                first_spike(i) = find(xe(index_e(i), :) == 1, 1);
        else
                first_spike(i) = 1e4;
        end
end

M = [(1:length(index_e))', first_spike];
sortM = sortrows(M, 2);
sort_ind_e = sortM(:, 1);
plot_spikes_e = xe(index_e(sort_ind_e(1:samples:length(index_e))), 1:206);
plot_first_spikes_e = zeros(size(plot_spikes_e));
plot_spikes_e_2 = zeros(size(plot_spikes_e));
for i = 1:size(plot_spikes_e, 1)
        if first_spike(sort_ind_e(1+(i-1)*samples)) < size(plot_spikes_e, 2) - 6
                plot_first_spikes_e(i, first_spike(sort_ind_e(1+(i-1)*samples))+(0:6)) = 1;
        end
        for j = 1:size(plot_spikes_e, 2)
                if plot_spikes_e(i, j) == 1
                        plot_spikes_e_2(i, j+(0:6)) = 1;
                end
        end   
end
plot_spikes_e = plot_spikes_e_2;
plot_other_spikes_e = plot_spikes_e - plot_first_spikes_e;
color_channels_e = ones([size(plot_spikes_e), 3]);
color_channels_e(:, :, 1) = color_channels_e(:, :, 1) - plot_first_spikes_e;
color_channels_e(:, :, 2) = color_channels_e(:, :, 2) - plot_first_spikes_e;
color_channels_e(:, :, 3) = color_channels_e(:, :, 3) - 0.3*plot_first_spikes_e;
color_channels_e(:, :, 1) = color_channels_e(:, :, 1) - plot_other_spikes_e/3;
color_channels_e(:, :, 2) = color_channels_e(:, :, 2) - plot_other_spikes_e/3;
color_channels_e(:, :, 3) = color_channels_e(:, :, 3) - 0.1*plot_other_spikes_e;

figure;
set(gcf, 'color', [1 1 1])
imagesc((1:size(color_channels_e, 2))*0.1-0.3, 1:size(color_channels_e, 1), color_channels_e)
xlabel('Time (ms)')
ylabel('#Neurons')
ax = gca;
a = get(ax,'Label');
set(gca,'Label', a,'fontsize',15);
savefig('sequence/e_sequence_visualize.fig');
saveas(gcf, 'sequence/e_sequence_visualize.pdf');

first_spike = zeros(length(index_i), 1);
for i = 1:length(index_i)
        if sum(xi(index_i(i), :)) > 0
                first_spike(i) = find(xi(index_i(i), :) == 1, 1);
        else
                first_spike(i) = 1e4;
        end
end

M = [(1:length(index_i))', first_spike];
sortM = sortrows(M, 2);
sort_ind_i = sortM(:, 1);
plot_spikes_i = xi(index_i(sort_ind_i(1:samples:length(index_i))), 1:206);
plot_first_spikes_i = zeros(size(plot_spikes_i));
plot_spikes_i_2 = zeros(size(plot_spikes_i));
for i = 1:size(plot_spikes_i, 1)
        if first_spike(sort_ind_i(1+(i-1)*samples)) < size(plot_spikes_i, 2) - 6 
                plot_first_spikes_i(i, first_spike(sort_ind_i(1+(i-1)*samples))+(0:6)) = 1;
        end
        for j = 1:size(plot_spikes_i, 2)
                if plot_spikes_i(i, j) == 1
                        plot_spikes_i_2(i, j+(0:6)) = 1;
                end
        end
end
plot_spikes_i = plot_spikes_i_2;
plot_other_spikes_i = plot_spikes_i - plot_first_spikes_i;
color_channels_i = ones([size(plot_spikes_i), 3]);
color_channels_i(:, :, 3) = color_channels_i(:, :, 3) - plot_first_spikes_i;
color_channels_i(:, :, 2) = color_channels_i(:, :, 2) - plot_first_spikes_i;
color_channels_i(:, :, 1) = color_channels_i(:, :, 1) - 0.3*plot_first_spikes_i;
color_channels_i(:, :, 3) = color_channels_i(:, :, 3) - plot_other_spikes_i/3;
color_channels_i(:, :, 2) = color_channels_i(:, :, 2) - plot_other_spikes_i/3;
color_channels_i(:, :, 1) = color_channels_i(:, :, 1) - 0.1*plot_other_spikes_i;

figure;
set(gcf, 'color', [1 1 1])
imagesc((1:size(color_channels_i, 2))*0.1-0.3, 1:size(color_channels_i, 1), color_channels_i)
xlabel('Time (ms)')
ylabel('#Neurons')
ax = gca;
a = get(ax,'Label');
set(gca,'Label', a,'fontsize',15);
savefig('sequence/i_sequence_visualize.fig');
saveas(gcf, 'sequence/i_sequence_visualize.pdf');

%{
color_channels_ei = ones([max(size(plot_spikes_e), size(plot_spikes_i)), 3]);
color_channels_ei(:, :, 3) = color_channels_ei(:, :, 3) - plot_first_spikes_i - 0.3*plot_first_spikes_e;
color_channels_ei(:, :, 2) = color_channels_ei(:, :, 2) - plot_first_spikes_i - plot_first_spikes_e;
color_channels_ei(:, :, 1) = color_channels_ei(:, :, 1) - 0.3*plot_first_spikes_i - plot_first_spikes_e;
color_channels_ei(:, :, 3) = color_channels_ei(:, :, 3) - plot_other_spikes_i/3 - 0.1* plot_other_spikes_e;
color_channels_ei(:, :, 2) = color_channels_ei(:, :, 2) - plot_other_spikes_i/3 - plot_other_spikes_e/3;
color_channels_ei(:, :, 1) = color_channels_ei(:, :, 1) - 0.1*plot_other_spikes_i - plot_other_spikes_i/3;

figure;
set(gcf, 'color', [1 1 1])
imagesc((1:size(color_channels_ei, 2))*0.1-0.3, 1:size(color_channels_ei, 1), color_channels_ei)
xlabel('Time (ms)')
ylabel('#Neurons')
ax = gca;
a = get(ax,'Label');
set(gca,'Label', a,'fontsize',15);
savefig('sequence/ei_sequence_visualize.fig');
saveas(gcf, 'sequence/ei_sequence_visualize.pdf');
%}
