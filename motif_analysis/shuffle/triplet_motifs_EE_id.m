%record the id of each motif

net_num = 10;

network_folder = '../../basic_network/2/';

V0 = -0.063;

ckpt_num = 2001;

th0 = 2.5;

motif_types = 5;

neurons = 1200;

num_pairs = neurons*(neurons - 1)/2;

num_triplets = neurons*(neurons-1)*(neurons-2)/6;

p0 = zeros(ckpt_num, net_num); %strong connections

n1 = zeros(ckpt_num, net_num);  %linear chain
n2 = zeros(ckpt_num, net_num);  %fan-out
n3 = zeros(ckpt_num, net_num);  %fan-in
n4 = zeros(ckpt_num, net_num);  %fan-in/out
n5 = zeros(ckpt_num, net_num);  %circle

th_e_all = cell(net_num, 1);

motif_id = cell(net_num, 1);

for net = 1:net_num
        motif_id{net} = cell(motif_types, ckpt_num);
        for j = 1:motif_types
                for ckpt = 1:ckpt_num
                        motif_id{net}{j, ckpt} = [];
                end
        end
end

tau_ee = 1719;  %from turnover rate
p_swap = 1 - exp(-1/tau_ee);

parfor net = 1:net_num
        
        network_dir = [network_folder, num2str(net), '/'];
        
        EE_seg0 = h5read_sparsematrix([network_dir, '/init_condition_w.h5'], '/EE');
        ET0 = h5read([network_dir, '/init_condition_w.h5'], '/ET');
        EE_seg0.nzval = EE_seg0.nzval ./ (ET0(EE_seg0.rowval) - V0)/1e3; %transform to effw
        
        EE_seg = EE_seg0;

        for ckpt = 1:ckpt_num
                
                disp(['Net ', num2str(net), ', checkpoint ', num2str(ckpt)]);
                
                effw = retrieve_sp(EE_seg);

                backbone = (effw > th0);

                p0(ckpt, net) = sum(sum(backbone))/num_pairs/2;
                
                for i = 1:neurons
                        for j = 1:neurons
                                if backbone(j, i) == 1 && j ~= i
                                        for k = 1:neurons
                                                if backbone(k, j) == 1 && k ~= i && k~= j
                                                        if backbone(k, i) == 0 && backbone(i, j) == 0 && ...
                                                                        backbone(i, k) == 0 && backbone(j, k) == 0
                                                                n1(ckpt, net) = n1(ckpt, net)+1;
                                                                motif_id{net}{1, ckpt} = [motif_id{net}{1, ckpt}; [i, j, k]];
                                                        end
                                                end
                                                
                                                if backbone(k, i) == 1 && k ~= i && k > j
                                                        if backbone(i, j) == 0 && backbone(i, k) == 0 && ...
                                                                        backbone(j, k) == 0 && backbone(k, j) == 0
                                                                n2(ckpt, net) = n2(ckpt, net)+1;
                                                                motif_id{net}{2, ckpt} = [motif_id{net}{2, ckpt}; [i, j, k]];
                                                        end
                                                end
                                                
                                                if backbone(k, i) == 1 && backbone(k, j) == 1 && k~= i && k ~= j
                                                        if backbone(i, j) == 0 && backbone(i, k) == 0 && backbone(j, k) == 0
                                                                n4(ckpt, net) = n4(ckpt, net)+1;
                                                                motif_id{net}{4, ckpt} = [motif_id{net}{4, ckpt}; [i, j, k]];
                                                        end
                                                end
                                                
                                                if backbone(i, k) == 1 && backbone(k, j) == 1 && k < i && k < j
                                                        if backbone(i, j) == 0 && backbone(k, i) == 0 && backbone(j, k) == 0
                                                                n5(ckpt, net) = n5(ckpt, net)+1;
                                                                motif_id{net}{5, ckpt} = [motif_id{net}{5, ckpt}; [i, j, k]];
                                                        end
                                                end
                                                
                                        end
                                end
                                
                                if backbone(i, j) == 1 && j ~= i
                                        for k = 1:neurons
                                                if backbone(i, k) == 1 && k ~= i && k > j
                                                        if backbone(j, i) == 0 && backbone(k, i) == 0 && ...
                                                                        backbone(j, k) == 0 && backbone(k, j) == 0
                                                                n3(ckpt, net) = n3(ckpt, net)+1;
                                                                motif_id{net}{3, ckpt} = [motif_id{net}{3, ckpt}; [i, j, k]];
                                                        end
                                                end
                                        end
                                end
           
                        end

                end
                
                strong_idx = find(EE_seg.nzval > th0);
                weak_idx = find(EE_seg.nzval <= th0);
                swap_idx = (rand(length(strong_idx), 1) < p_swap);
                swap_strong_idx = strong_idx(swap_idx == 1);
                swap_weak_idx = weak_idx(randperm(length(weak_idx), sum(swap_idx)));
                temp = EE_seg.nzval(swap_weak_idx);
                EE_seg.nzval(swap_weak_idx) = EE_seg.nzval(swap_strong_idx);
                EE_seg.nzval(swap_strong_idx) = temp;
                
        end
end

n = cat(3, n1, n2, n3, n4, n5);

E_num = zeros(ckpt_num, net_num, motif_types);
cdf_all = zeros(ckpt_num, net_num, motif_types);

mu = zeros(ckpt_num, motif_types);
sigma = zeros(ckpt_num, motif_types);
p_value = zeros(ckpt_num, motif_types);

for ckpt = 1:ckpt_num
        
        p_all = zeros(motif_types, net_num);
        
        p_all(1, :) = 6*p0(ckpt, :).^2 .*(1-p0(ckpt, :)).^4;
        p_all(2, :) = 3*p0(ckpt, :).^2 .*(1-p0(ckpt, :)).^4;
        p_all(3, :) = 3*p0(ckpt, :).^2 .*(1-p0(ckpt, :)).^4;
        p_all(4, :) = 6*p0(ckpt, :).^3 .*(1-p0(ckpt, :)).^3;
        p_all(5, :) = 2*p0(ckpt, :).^3 .*(1-p0(ckpt, :)).^3;
        
        for motif = 1:motif_types
                E_num(ckpt, :, motif) =  p_all(motif, :) * num_triplets;
                mu(ckpt, motif) = sum(E_num(ckpt, :, motif));
                sigma(ckpt, motif) = sqrt(sum(E_num(ckpt, :, motif) .* (1-p_all(motif, :))));
                p_value(ckpt, motif) = cdf('normal', (mu(ckpt, motif) - sum(n(ckpt, :, motif)))/sigma(ckpt, motif), 0, 1);
        end
        
        for net = 1:net_num
                for motif = 1:motif_types
                        cdf_all(ckpt, net, motif) = cdf('Binomial', n(ckpt, net, motif), num_triplets, p_all(motif, net));
                end
        end

end

save('triplet_motifs_EE_id.mat', 'n', 'E_num', 'cdf_all', 'mu', 'sigma', 'p_value', 'p0', 'motif_id', '-v7.3');

