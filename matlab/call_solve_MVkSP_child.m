function result = call_solve_MVkSP_child(A, k, F0, F1, gub, x_gub, init_B_TRI, init_l_TRI, params)

    %disp(params)
    
    %disp("F0")
    %disp(F0)
    %disp("F1")
    %disp(F1)

    %disp(size(init_B_TRI))
    %disp(size(init_l_TRI))

    result = solve_MVkSP_sdpnalplus(A, k, F0, F1, gub, x_gub, init_B_TRI, init_l_TRI, params);

    if result.best_gap <= params.bb_tol
        result.br_idx = -1;
    else
        result.br_idx = find_branching_index(result.best_x_lb);
    end

    %disp(result)


end