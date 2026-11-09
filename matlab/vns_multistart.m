function [f_best, x_best] = vns_multistart(Q, c, k_card, n_times)

    x_best = [];
    f_best = inf;
    n = size(Q, 1);

    VNS.iter_max = 2;
    VNS.kmin     = 2;
    VNS.kmax     = floor(n/2);
    VNS.kstep    = 1;
    VNS.max_time_sec = 30;
    VNS.verbose  = false;
    
    for i=1:n_times

        x0 = randn(n,1);
        VNS.seed = i;
        [f_vns, x_vns] = vns_binary_k_qp(Q, c, x0, k_card, VNS);
         if f_vns < f_best
             disp(f_vns)
             f_best = f_vns;
             x_best = x_vns;
         end
    
    end

end