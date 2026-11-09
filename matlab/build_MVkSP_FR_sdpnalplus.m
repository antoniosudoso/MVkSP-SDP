function [blk, At, b, C, L, U] = build_MVkSP_FR_sdpnalplus(Q, c, r)
    % Build SDPNAL+ data for
    %
    %   min_Z  c' diag(Z) + <Q,Z>
    %   s.t.   diag(Z) = (1/r) Z*1
    %          tr(Z) = r
    %          Z >= 0
    %          Z >= 0 elementwise
    %
    % with Z in S^|U|.

    [~, n] = size(Q);

    % One PSD block of size |U|
    blk = {'s', n};

    C = cell(1);
    C{1} = spdiags(c, 0, n, n) + Q;

    % n constraints from diag(Z) = (1/r) Z*1
    % 1 constraint from 1'Z1 = r^2
    Acell = cell(1, n + 1);
    b = zeros(n + 1, 1);

    for i = 1:n
        % Build A_i such that
        % <A_i, Z> = Z_ii - (1/r) * sum_j Z_ij

        rows = [i, repmat(i, 1, n-1), setdiff(1:n, i)];
        cols = [i, setdiff(1:n, i), repmat(i, 1, n-1)];
        vals = [(1 - 1/r), repmat(-1/(2*r), 1, 2*(n-1))];

        Acell{i} = sparse(rows, cols, vals, n, n);
        b(i) = 0;
    end

    % Last constraint: 1' Z 1 = r^2
    %Acell{n+1} = sparse(ones(n, n));
    %b(n+1) = r^2;
    Acell{n+1} = speye(n);
    b(n+1) = r;

    % Convert to SDPNAL+ svec format
    At = svec(blk, Acell, 1);

    L = cell(1);
    L{1} = 0;
    U = [];

end