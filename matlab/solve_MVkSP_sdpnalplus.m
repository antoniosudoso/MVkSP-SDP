function result = solve_MVkSP_sdpnalplus(A, k, F0, F1, gub, x_gub, init_B_TRI, init_l_TRI, params)

    tStartCP = tic;

    maxNumCompThreads(params.n_threads);

    [~, original_n] = size(A);
    [U_idx, Q, c, r, const] = get_QP_data(A, k, F0, F1);
    [~, n] = size(Q);

    %keyboard

    VNS.iter_max = 2;
    VNS.kmin     = 2;
    VNS.kmax     = floor(n/2);
    VNS.kstep    = 2;
    VNS.max_time_sec = 30;
    VNS.verbose  = false;

    result = struct();
    result.n = n;
    
    result.ineq_list = [];
    result.timeSDP_list = [];
    result.lb_list = [];
    result.ub_list = [];
    result.gap_list = [];
    
    result.best_lb = -inf;
    result.best_ub = gub;
    result.best_gap = inf;
    result.best_x_ub = x_gub;
    result.best_x_lb = [];
    result.cp_flag = -1;
    result.cp_iter = 0;
    result.time = 0;

    if params.sdp_fr
        [blk, At, b, C, L, U] = build_MVkSP_FR_sdpnalplus(Q, c, r);
        max_eig = r;
    else
        [blk, At, b, C, L, U] = build_MVkSP_sdpnalplus(Q, c, r);
        max_eig = r + 1;
    end


    if params.sdp_fr
        [init_B_TRI, init_l_TRI] = reduce_triangle(F0, F1, init_B_TRI, init_l_TRI, n, original_n);
    else
        [init_B_TRI, init_l_TRI] = reduce_triangle(F0 + 1, F1 + 1, init_B_TRI, init_l_TRI, n + 1, original_n + 1);
    end

    result.best_Bcell = init_B_TRI;
    result.best_l = init_l_TRI;

    %keyboard

    OPTIONS.printlevel = params.sdp_verbose;
    OPTIONS.tol = params.sdp_tol;
    OPTIONS.stopoption = 0;
    OPTIONS.maxiter = 100000;
    OPTIONS.maxtime = 3600;

    Bt = svec(blk, init_B_TRI, 1);
    l = init_l_TRI;

    % ----- Solve -----
    tStart = tic;
    [~, Y, ~, y, S, Z, y2, ~, ~, ~] = sdpnalplus(blk, At, C, b, L, U, Bt, l, [], OPTIONS);
    timeSDP = toc(tStart);
    [~, lb_safe] = safe_bound_error(blk, At, C, b, y, Z, Bt, y2, l, max_eig);
    lb_safe = lb_safe + const;
    [~, lb_lp, timeLP] = safe_bound_lp(blk, S, At, Bt, C, b, l, params.n_threads);
    lb_lp = lb_lp + const;
    lb = max(lb_safe, lb_lp);
    
    result.best_lb = lb;
    result.lb_list = [result.lb_list; lb];
    result.timeSDP_list = [result.timeSDP_list; timeSDP];
    result.ineq_list = [result.ineq_list; size(init_B_TRI, 2)];

    % Extract solution to initialize the heuristic
    Y_SDP = Y{1};

    if params.sdp_fr
        x_SDP = diag(Y_SDP);
    else
        x_SDP = Y_SDP(2:n+1, 1);
    end
    % Recover full x
    x_lb = zeros(original_n, 1);
    x_lb(F1) = 1;
    x_lb(U_idx) = x_SDP;
    result.best_x_lb = x_lb;

    [f_vns, x_vns] = vns_binary_k_qp(Q, c, x_SDP, r, VNS);     
    f_vns = f_vns + const;
    if f_vns < result.best_ub
        result.best_ub = f_vns;
        % Recover full x
        x_star = zeros(original_n, 1);
        x_star(F1) = 1;
        x_star(U_idx) = x_vns;
        result.best_x_ub = x_star;
    end

    result.ub_list = [result.ub_list; f_vns];
    gap = (result.best_ub - lb) / max(1, abs(result.best_ub));
    result.gap_list = [result.gap_list; gap];
    result.best_gap = gap;

    fprintf('\n\t ************************************\n')
    fprintf('\t Iteration             = %d \n',    result.cp_iter);
    fprintf('\t Inequalities          = %d \n',    size(init_B_TRI, 2)); 
    fprintf('\t Valid LB              = %.4f \n',  lb);
    fprintf('\t LB error bound        = %.4f \n',  lb_safe);
    fprintf('\t LB LP                 = %.4f \n',  lb_lp);
    fprintf('\t Time SDP (s)          = %.4f \n',  timeSDP);
    fprintf('\t Time LP (s)           = %.4f \n',  timeLP);
    fprintf('\t UB                    = %.4f \n',  result.best_ub);
    fprintf('\t Relative gap (%%)      = %.4f \n',   result.best_gap*100);
    fprintf('\t ************************************\n\n')

    if result.best_gap <= params.bb_tol
        result.cp_flag = 0;
        result.time = toc(tStartCP);
        disp(' ');
        return
    end

    B_TRI = init_B_TRI;
    l_TRI = init_l_TRI;

    for i=1:params.cp_maxiter
        
        if ~isempty(B_TRI) && ~isempty(l_TRI)
            Bt = svec(blk, B_TRI, 1);
            Yvec = svec(blk, Y{1}, 1);
            active_id = abs(Bt{1}' * Yvec - l_TRI) <= params.cp_activeineq;
            B_TRI = B_TRI(active_id);
            l_TRI = l_TRI(active_id);
        end
        
        if params.sdp_fr
            [B_triangle, l_triangle, ~] = separate_RLT_triangle_FR(Y{1}, params.cp_epsineq, params.cp_maxineq, params.cp_newineq);
        else
            [B_triangle, l_triangle, ~] = separate_RLT_triangle(Y{1}, params.cp_epsineq, params.cp_maxineq, params.cp_newineq);
        end
        n_triangle = size(B_triangle, 2);
        if n_triangle < n
            result.cp_flag = 1;
            disp(' ');
            break;
        end

        B_TRI = [B_TRI, B_triangle];
        l_TRI = [l_TRI; l_triangle];
        Bt = svec(blk, B_TRI, 1);
        
        % solve SDP
        tStart = tic;
        [~, Y, ~, y, S, Z, y2, ~, ~, ~] = sdpnalplus(blk, At, C, b, L, U, Bt, l_TRI, [], OPTIONS);
        timeSDP = toc(tStart);
        [~, lb_safe] = safe_bound_error(blk, At, C, b, y, Z, Bt, y2, l_TRI, max_eig);
        lb_safe = lb_safe + const;
        [~, lb_lp, timeLP] = safe_bound_lp(blk, S, At, Bt, C, b, l_TRI, params.n_threads);
        lb_lp = lb_lp + const;
        lb = max(lb_safe, lb_lp);

        result.ineq_list = [result.ineq_list; size(B_TRI, 2)];
        result.lb_list = [result.lb_list; lb];
        result.timeSDP_list = [result.timeSDP_list; timeSDP];

        % Extract solution to initialize the heuristic
        Y_SDP = Y{1};

        if params.sdp_fr
            x_SDP = diag(Y_SDP);
        else
            x_SDP = Y_SDP(2:n+1, 1);
        end

        [f_vns, x_vns] = vns_binary_k_qp(Q, c, x_SDP, r, VNS);   
        f_vns = f_vns + const;
        if f_vns < result.best_ub
            result.best_ub = f_vns;
            % Recover full x
            x_star = zeros(original_n, 1);
            x_star(F1) = 1;
            x_star(U_idx) = x_vns;
            result.best_x_ub = x_star;
        end
    
        result.ub_list = [result.ub_list; f_vns];
        gap = (result.best_ub - lb) / max(1, abs(result.best_ub));
        result.gap_list = [result.gap_list; gap];

        result.cp_iter = result.cp_iter + 1;

        if lb > result.best_lb
            result.best_lb = lb;
            result.best_gap = gap;
            % Recover full x
            x_lb = zeros(original_n, 1);
            x_lb(F1) = 1;
            x_lb(U_idx) = x_SDP;
            result.best_x_lb = x_lb;
            result.best_Bcell = B_TRI;
            result.best_l = l_TRI;
        else
            result.cp_flag = 2;
            disp(' ');
            break
        end
    
        fprintf('\n\t ************************************\n')
        fprintf('\t Iteration             = %d \n',    result.cp_iter);
        fprintf('\t Inequalities          = %d \n',    size(B_TRI, 2)); 
        fprintf('\t Valid LB              = %.4f \n',  lb);
        fprintf('\t LB error bound        = %.4f \n',  lb_safe);
        fprintf('\t LB LP                 = %.4f \n',  lb_lp);
        fprintf('\t Time SDP (s)          = %.4f \n',  timeSDP);
        fprintf('\t Time LP (s)           = %.4f \n',  timeLP);
        fprintf('\t UB                    = %.4f \n',  result.best_ub);
        fprintf('\t Relative gap (%%)      = %.4f \n',   result.best_gap*100);
        fprintf('\t ************************************\n\n')

        if result.best_gap <= params.bb_tol
            result.cp_flag = 0;
            disp(' ');
            break
        end

    end

    if params.sdp_fr
        [or_B_TRI, or_l_TRI] = restore_triangle(F0, F1, result.best_Bcell, result.best_l, original_n);
    else
        [or_B_TRI, or_l_TRI] = restore_triangle(F0 + 1, F1 + 1, result.best_Bcell, result.best_l, original_n + 1);
    end

    result.best_Bcell = or_B_TRI;
    result.best_l = or_l_TRI;

    result.time = toc(tStartCP);

end