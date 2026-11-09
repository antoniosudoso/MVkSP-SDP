function [U, Q, c, r, const] = get_QP_data(A, k, F0, F1)

    [d, n] = size(A);

    % Ensure column vectors, unique sorted indices
    F0 = unique(F0(:));
    F1 = unique(F1(:));

    % Basic validity checks
    if any(F0 < 1) || any(F0 > n) || any(F1 < 1) || any(F1 > n)
        error('F0 and F1 must contain indices between 1 and n.');
    end
    if ~isempty(intersect(F0,F1))
        error('F0 and F1 must be disjoint.');
    end

    s = numel(F1);
    r = k - s;

    % Remaining unfixed indices
    U = setdiff((1:n)', union(F0,F1));

    % Infeasible node checks
    if r < 0
        error('Infeasible node: |F1| > k.');
    end
    if r > numel(U)
        error('Infeasible node: not enough unfixed variables to reach cardinality k.');
    end

    % Precompute full norms
    c = sum(A.^2,1)';   % c_i = ||p_i||^2

    % Fixed-to-1 contribution q = sum_{i in F1} p_i
    if isempty(F1)
        q = zeros(d,1);
        const = 0;
    else
        q = sum(A(:,F1), 2);
        const = (1/k) * sum(c(F1)) - (1/k^2) * (q' * q);
    end

    % Submatrix for unfixed indices
    A_U = A(:,U);
    c_U = c(U);
    G_U = A_U' * A_U;

    % Node objective: const + c' * y + (y' * Q * y);

    % Linear term c
    c = (1/k) * c_U - (2/k^2) * (A_U' * q);
    % Quadratic term Q
    Q = - (1/(k^2)) * G_U;

end