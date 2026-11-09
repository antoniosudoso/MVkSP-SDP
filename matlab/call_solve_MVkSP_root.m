function result = call_solve_MVkSP_root(A, k, params)

    disp(params)

    Q = -(1/(k^2))*(A'*A);
    c = (1/k)*sum(A.^2,1)';

    [f_vns, x_vns] = vns_multistart(Q, c, k, 10);

    F0 = [];
    F1 = [];
    init_B_TRI = {};
    init_l_TRI = [];
    gub = f_vns;
    x_gub = x_vns;

    result = solve_MVkSP_sdpnalplus(A, k, F0, F1, gub, x_gub, init_B_TRI, init_l_TRI, params);

    if result.best_gap <= params.bb_tol
        result.br_idx = -1;
    else
        result.br_idx = find_branching_index(result.best_x_lb);
    end

    disp(result)


end