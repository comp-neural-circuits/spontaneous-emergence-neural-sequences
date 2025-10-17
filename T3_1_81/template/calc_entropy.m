%run after delay_and_jitter.m

delay_all_neurons = [delay_e_all; delay_i_all];
delay_all_neurons(sum(abs(delay_all_neurons), 2) == ...
        min(sum(abs(delay_all_neurons), 2)), :) = [];
%remove the source

norm_entropy_data = normalized_entropy(delay_all_neurons);

shuffle_num = 10;

norm_entropy_shuffle = zeros(shuffle_num, length(norm_entropy_data));

for shuffle = 1:shuffle_num
        delay_shuffle = zeros(size(delay_all_neurons));
        for trial = 1:size(delay_all_neurons, 2)
                delay_shuffle(:, trial) = delay_all_neurons(randperm(size(delay_all_neurons, 1)), trial);
        end
        norm_entropy_shuffle(shuffle, :) = normalized_entropy(delay_shuffle);
end

save('sequence/entropy.mat', 'norm_entropy_data', 'norm_entropy_shuffle')

function entropy_norm = normalized_entropy(delay_all_neurons)

rank_num = size(delay_all_neurons, 1);
trials = size(delay_all_neurons, 2);
neuron_rank = zeros(size(delay_all_neurons));

for trial = 1:trials
        for neuron = 1:rank_num
                if delay_all_neurons(neuron, trial) == -1
                        neuron_rank(neuron, trial) = -1;
                else
                        neuron_rank(neuron, trial) = sum(delay_all_neurons(:, trial) <= ...
                        delay_all_neurons(neuron, trial) & delay_all_neurons(:, trial) > 0);
                end
        end
end
                
prob_all = zeros(rank_num, rank_num);

for k = 1:rank_num
        for neuron = 1:rank_num
                prob_all(neuron, k) = mean(neuron_rank(neuron, :) == k);
        end
end

prob_all = prob_all + ones(rank_num, 1)*(1 - sum(prob_all))/(rank_num);

entropy_k = zeros(1, rank_num);

for k = 1:rank_num
        for neuron = 1:rank_num
                if prob_all(neuron, k) > 0
                        entropy_k(k) = entropy_k(k) - prob_all(neuron, k)*log2(prob_all(neuron, k));
                end
        end
end

entropy_norm = entropy_k/log2(rank_num);

end

