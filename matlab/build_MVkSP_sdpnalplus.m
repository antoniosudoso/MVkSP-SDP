function [blk, At, b, C, L, U] = build_MVkSP_sdpnalplus(Q, c, r)
    % Build SDPNAL+ data for the DNN relaxation
    %
    % Variable:
    %   Y = [ 1   x';
    %         x    X ]  in S^{|U|+1}
    %
    % Objective:
    %   min c'x + <Q,X>
    %
    % Constraints:
    %   1'x = r
    %   diag(X) = x
    %   X*1 = r x
    %   Y >= 0 elementwise
    %   x <= 1 (upper bound only on Y(1,i+1), Y(i+1,1))
    %   Y(1,1) = 1

    [~, n] = size(Q);

    % ----- SDPNAL+ blocks -----
    msize = n + 1;      % Y is (n+1)x(n+1)
    blk = cell(1,2);
    blk{1,1} = 's';
    blk{1,2} = msize;

    % ----- Objective: <C, Y> -----
    % x_i = Y(1,i+1)
    % X_ij = Y(i+1,j+1)
    Cmat = sparse(msize, msize);
    Cmat(2:end, 2:end) = Q;

    % Put (1/k)*c'x on off-diagonal (1,i+1) and (i+1,1)
    Cmat(1, 2:end) = (1/2) * c';
    Cmat(2:end, 1) = (1/2) * c;

    C = {Cmat};

    % ----- Bounds L <= Y <= U -----
    % DNN nonnegativity: Y >= 0 elementwise
    L = {0};
    U = [];

    % ----- Equality constraints <A_t, Y> = b_t -----
    Acell = {};
    b = [];

    % (0) Y(1,1) = 1
    A00 = sparse(msize, msize);
    A00(1,1) = 1;
    Acell{end+1} = A00;
    b(end+1,1) = 1;

    % (1) sum(x) = r, where x_i = Y(1,i+1)
    A1 = zeros(msize, msize);
    for i = 1:n
        A1(1, i+1) = A1(1, i+1) + 0.5;
        A1(i+1, 1) = A1(i+1, 1) + 0.5;
    end
    Acell{end+1} = sparse(A1);
    b(end+1,1) = r;

    % (2) diag(X) = x -> Y(i+1,i+1) - Y(1,i+1) = 0
    for i = 1:n
        Ai = zeros(msize, msize);
        Ai(i+1, i+1) = 1;
        Ai(1, i+1)   = Ai(1, i+1) - 0.5;
        Ai(i+1, 1)   = Ai(i+1, 1) - 0.5;
        Acell{end+1} = sparse(Ai);
        b(end+1,1) = 0;
    end

    % (3) X*1 = k x
    % For each i:
    %   sum_{j=1..n} X(i,j) - k*x_i = 0
    % In Y indexing:
    %   sum_{j=1..n} Y(i+1,j+1) - k*Y(1,i+1) = 0
    for i = 1:n
        Ai = zeros(msize, msize);

        % Add sum_j X(i,j):
        % diagonal term X(i,i) has coefficient 1 -> Y(i+1,i+1) += 1
        Ai(i+1, i+1) = Ai(i+1, i+1) + 1;

        % off-diagonals X(i,j), j!=i: coefficient 1 on Y(i+1,j+1)
        % represent symmetric inner product by placing 0.5 on both (i+1,j+1) and (j+1,i+1)
        for j = 1:n
            if j == i, continue; end
            Ai(i+1, j+1) = Ai(i+1, j+1) + 0.5;
            Ai(j+1, i+1) = Ai(j+1, i+1) + 0.5;
        end

        % Subtract k*x_i: x_i = Y(1,i+1)
        Ai(1, i+1)   = Ai(1, i+1) - (r/2);
        Ai(i+1, 1)   = Ai(i+1, 1) - (r/2);

        Acell{end+1} = sparse(Ai);
        b(end+1,1) = 0;
    end

    % Convert to At using SDPNAL+ svec format
    At = svec(blk, Acell, 1);

end