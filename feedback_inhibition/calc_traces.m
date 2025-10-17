control_file = 'control/testing_control.mat';
test_dir = 'source/';

load(control_file)
rcd0 = rcd;
avg_V0 = mean(rcd0.Ve(:, 1501:2000), 2);
files = dir([test_dir, 'testing_inject_point_*.mat']);

Ve = zeros(length(files), size(rcd.xe, 2));
avg_V = zeros(length(files), 1);

for i = 1:length(files)
        inject_point = str2double(files(i).name(22:end-4));
        load([test_dir, files(i).name]);
        Ve(i, :) = rcd.Ve(inject_point, :);
        avg_V(i) = mean(rcd.Ve(inject_point, 1501:2000));
end

Ve0_avg = mean(rcd0.Ve);
Ve0_std = std(rcd0.Ve);

save('isolated_and_source_neurons.mat', 'avg_V0', 'Ve0_avg', 'Ve0_std', 'Ve', 'avg_V');