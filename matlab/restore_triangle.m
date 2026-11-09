function [or_B_TRI, or_l_TRI] = restore_triangle(F0, F1, shr_B_TRI, shr_l_TRI, original_n)
   

    n_ineq = size(shr_B_TRI, 2);
    or_B_TRI = cell(1, n_ineq);
    or_l_TRI = zeros(n_ineq, 1);
        
    counter = 1;

    x_id = setdiff(1:original_n, union(F0, F1), 'stable');

    for c=1:n_ineq
        
        [id_i, id_j, v] = find(shr_B_TRI{c});
        or_B_TRI{counter} = sparse(x_id(id_i), x_id(id_j), v, original_n, original_n);
        or_l_TRI(counter) = shr_l_TRI(c);
        counter = counter + 1;
        
    end
    
end