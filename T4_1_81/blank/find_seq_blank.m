T = 43;
dt = 1e-4;
frames = round(T/dt);

bin_width = 0.1;
bin_num = round(T/bin_width);

Ne = 1200;
Ni = 240;

for sim = 1:10

disp(sim)
folder = num2str(sim);

re = zeros(Ne, bin_num);
ri = zeros(Ni, bin_num);

parfor bin = 1:bin_num
        xe_seg = h5read_sparsematrix([folder, '/all_x.h5'], [num2str(bin), '/xe']);
        xi_seg = h5read_sparsematrix([folder, '/all_x.h5'], [num2str(bin), '/xi']);
        xe = retrieve_sp(xe_seg);
        xi = retrieve_sp(xi_seg);
        re(:, bin) = sum(xe(:, 1:end-1), 2)/bin_width;
        ri(:, bin) = sum(xi(:, 1:end-1), 2)/bin_width;
end

stim_onset = 3.1;
stim_period = 0.4;
stim_num = 100;

bin_baseline = round(stim_onset/bin_width):round(stim_period/bin_width):...
        round((stim_onset+(stim_num-1)*stim_period)/bin_width);
bin_effect = [bin_baseline+1; bin_baseline+2; bin_baseline+3];

re_baseline = zeros(Ne, stim_num);
ri_baseline = zeros(Ni, stim_num);
re_effect = zeros(Ne, stim_num);
ri_effect = zeros(Ni, stim_num);

for num = 1:stim_num
        re_baseline(:, num) = mean(re(:, bin_baseline(:, num)), 2);
        ri_baseline(:, num) = mean(ri(:, bin_baseline(:, num)), 2);
        re_effect(:, num) = mean(re(:, bin_effect(:, num)), 2);
        ri_effect(:, num) = mean(ri(:, bin_effect(:, num)), 2);
end

frm_e = mean(re_effect - re_baseline, 2);
frm_i = mean(ri_effect - ri_baseline, 2);

re_avg_baseline = mean(mean(re_baseline));
ri_avg_baseline = mean(mean(ri_baseline));

M = 1e3;    %number of null models
p = 0.95;

baseline = 0.1;
effect = 0.3;

fr_baseline_e = poissrnd(re_avg_baseline*baseline, M, stim_num)/baseline;
fr_effect_e = poissrnd(re_avg_baseline*effect, M, stim_num)/effect;
frm_null_e = mean(fr_effect_e - fr_baseline_e, 2);

fr_baseline_i = poissrnd(ri_avg_baseline*baseline, M, stim_num)/baseline;
fr_effect_i = poissrnd(ri_avg_baseline*effect, M, stim_num)/effect;
frm_null_i = mean(fr_effect_i - fr_baseline_i, 2);

thresh_e = prctile(frm_null_e, p*100);
thresh_i = prctile(frm_null_i, p*100);

follower_index_e = find(frm_e > thresh_e);
follower_index_i = find(frm_i > thresh_i);

max_frm_null_e = max(frm_null_e);
max_frm_null_i = max(frm_null_i);

save(['sequence_blank_', folder, '.mat'], 'follower_index_e', 'follower_index_i', 'frm_e', 'frm_i', 're_baseline', ...
        're_effect', 'ri_baseline', 'ri_effect', 'thresh_e', 'thresh_i', 'M', 'p', 'stim_onset', 'stim_period', 'stim_num', ...
        'T', 'baseline', 'effect', 'max_frm_null_e', 'max_frm_null_i');

end





