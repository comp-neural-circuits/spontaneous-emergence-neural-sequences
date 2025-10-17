%run after delay_and_jitter.m

%load('sequence/delay_and_jitter.mat')

if length(index_e) + length(index_i) > 1

        delay_all_neurons = [delay_e_all; delay_i_all];
        delay_all_neurons(sum(abs(delay_all_neurons), 2) == ...
                min(sum(abs(delay_all_neurons), 2)), :) = [];
        %remove the source

        med_delay = [med_delay_e; med_delay_i];
        med_delay(med_delay == min(med_delay)) = [];

        M = [(1:length(med_delay))', med_delay];
        sortM = sortrows(M, 2);
        sorted_neuron_idx = sortM(:, 1);

        r_data = calc_rank_corr(delay_all_neurons);

        shuffle_num = 10;

        r_shuffle = zeros(shuffle_num, length(r_data));

        for shuffle = 1:shuffle_num
                delay_shuffle = zeros(size(delay_all_neurons));
                for trial = 1:size(delay_all_neurons, 2)
                        delay_shuffle(:, trial) = delay_all_neurons(randperm(size(delay_all_neurons, 1)), trial);
                end
                r_shuffle(shuffle, :) = calc_rank_corr(delay_shuffle);
        end
        
else
        
        r_data = [];
        r_shuffle = [];
        sorted_neuron_idx = [];
        
end

save('sequence/rank_corr.mat', 'r_data', 'r_shuffle', 'sorted_neuron_idx');
                                

function r = calc_rank_corr(delay_all_neurons)

neurons = size(delay_all_neurons, 1);
trials = size(delay_all_neurons, 2);

O = zeros(neurons, neurons, trials);

for k = 1:trials
        for i = 1:neurons
                for j = 1:neurons
                        if  ((delay_all_neurons(i, k) > 0) && (delay_all_neurons(j, k) == -1)) || ...
                                ((delay_all_neurons(i, k) > 0) && (delay_all_neurons(j, k) > 0) && ...
                                (delay_all_neurons(i, k) < delay_all_neurons(j, k) - 2e-3))
                                  O(i, j, k) = 1;
                        end
                        if  ((delay_all_neurons(i, k) == -1) && (delay_all_neurons(j, k) > 0)) || ...
                                ((delay_all_neurons(i, k) > 0) && (delay_all_neurons(j, k) > 0) && ...
                                (delay_all_neurons(i, k) > delay_all_neurons(j, k) + 2e-3))
                                  O(i, j, k) = -1;
                        end
                end
        end
end

S = zeros(trials, trials, neurons);

parfor k = 1:trials
        disp(['k = ', num2str(k)])
        for l = 1:trials
                O_prod = O(:, :, k) .* O(:, :, l);
                for i = 1:neurons
                        S(k, l, i) = sum(O_prod(i, :))/(neurons - 1);
                end
        end
end

r = squeeze(mean(mean(S,  1), 2));

end
                                
                                
                                
                                
                                
                                