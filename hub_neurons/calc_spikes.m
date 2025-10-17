Ne = 1200;
Ni = 240;
trials = 10;
T = 6000;
spikes_e = cell(T, 1);
spikes_i = cell(T, 1);

for trial = 1:trials
        disp(trial)
        parfor i = 1:T
                disp(i)
                xe_seg = h5read_sparsematrix([num2str(trial), '/all_x.h5'], [num2str(i), '/xe']);
                xe = retrieve_sp(xe_seg);
                spikes_e{i, 1} = sum(xe, 2);
                xi_seg = h5read_sparsematrix([num2str(trial), '/all_x.h5'], [num2str(i), '/xi']);
                xi = retrieve_sp(xi_seg);
                spikes_i{i, 1} = sum(xi, 2);
        end
        save([num2str(trial), '/spikes_count.mat'], 'spikes_e', 'spikes_i')
end


