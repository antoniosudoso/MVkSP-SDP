function index = find_branching_index(x)

    n = length(x);
    pairs = sortrows([(1:n)', min(x, 1-x)], 2, 'descend'); % pairs (id_x, viol_bin)
    % disp(pairs)
    value = pairs(1, 2);
    if value <= 1e-6
        index = -1;
    else
        index = pairs(1, 1);
    end

end