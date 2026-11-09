function [best_cost, x_best] = vns_binary_k_qp(A, a, x0, k_card, opts)
% VNS_BINARY_K_QP
%   Variable Neighborhood Search for
%       min  x' * A * x + a' * x
%       s.t. x_i in {0,1}
%            sum(x) = k
%
% INPUT
%   A      : (n x n) symmetric matrix
%   a      : (n x 1) vector
%   x0     : (n x 1) initial vector (projected to {0,1} with sum = k_card)
%   k_card : required number of ones
%   opts   : struct with optional fields:
%              .iter_max        (default 2)
%              .kmin            (default 1)
%              .kmax            (default 100)
%              .kstep           (default 1)
%              .seed            (default 1)
%              .max_local_iters (default 1e4)
%              .max_time_sec    (default inf)
%              .verbose         (default false)
%
% OUTPUT
%   x_best    : best solution found (n x 1, entries in {0,1}, sum = k_card)
%   best_cost : objective value at x_best

    if nargin < 5
        opts = struct();
    end

    % ----- Options with defaults -----
    iter_max        = get_opt(opts, 'iter_max',        2);
    kmin            = get_opt(opts, 'kmin',            1);
    kmax            = get_opt(opts, 'kmax',            100);
    kstep           = get_opt(opts, 'kstep',           1);
    seed            = get_opt(opts, 'seed',            1);
    max_local_iters = get_opt(opts, 'max_local_iters', 1e4);
    max_time_sec    = get_opt(opts, 'max_time_sec',    60);
    verbose         = get_opt(opts, 'verbose',         false);

    % ----- Basic checks -----
    [n1, n2] = size(A);
    if n1 ~= n2
        error('A must be square.');
    end
    n = n1;

    if numel(a) ~= n
        error('Dimension mismatch: a must have length size(A,1).');
    end
    if numel(x0) ~= n
        error('Dimension mismatch: x0 must have length size(A,1).');
    end
    if ~isscalar(k_card) || ~isfinite(k_card) || k_card < 0 || k_card > n || floor(k_card) ~= k_card
        error('k_card must be an integer in [0, n].');
    end

    if ~isnumeric(A) || ~isreal(A) || any(~isfinite(A(:)))
        error('A must be a real finite numeric matrix.');
    end
    if ~isnumeric(a) || ~isreal(a) || any(~isfinite(a(:)))
        error('a must be a real finite numeric vector.');
    end
    if ~isnumeric(x0) || ~isreal(x0) || any(~isfinite(x0(:)))
        error('x0 must be a real finite numeric vector.');
    end

    if ~isscalar(seed) || ~isfinite(seed)
        error('seed must be a finite scalar.');
    end
    if ~isscalar(max_time_sec) || ~isreal(max_time_sec) || max_time_sec <= 0
        error('max_time_sec must be a positive finite scalar or inf.');
    end

    A = double(A);
    a = double(a(:));

    % ----- Initial feasible solution in {0,1} with sum = k_card -----
    x = project_to_binary_k(double(x0(:)), k_card);

    % ----- Initialize aggregates -----
    alphaA = A * x;      % A*x
    f1     = x' * A * x; % quadratic term
    f2     = a' * x;     % linear term

    x_best    = x;
    alphaA_b  = alphaA;
    f1_b      = f1;
    f2_b      = f2;
    best_cost = f1 + f2;

    % ----- RNG -----
    rng(seed);

    EPS = 1e-4;
    t0 = tic;

    % ===================== MAIN VNS LOOP =====================
    for it = 1:iter_max
        if toc(t0) > max_time_sec
            if verbose
                fprintf('Stopping: time limit reached in outer loop.\n');
            end
            break;
        end

        k = kmin;
        while k <= kmax
            if toc(t0) > max_time_sec
                if verbose
                    fprintf('Stopping: time limit reached in neighborhood loop.\n');
                end
                break;
            end

            % ---- Shake: perform k random feasible swaps ----
            [x, f1, f2, alphaA] = shake(x, f1, f2, alphaA, A, a, k);

            % ---- Local search: best improving swap ----
            [x, f1, f2, alphaA] = local_search_swap( ...
                x, f1, f2, alphaA, A, a, max_local_iters, max_time_sec, t0);

            cur_cost = f1 + f2;

            if cur_cost + EPS < best_cost
                x_best    = x;
                alphaA_b  = alphaA;
                f1_b      = f1;
                f2_b      = f2;
                best_cost = cur_cost;

                if verbose
                    fprintf('Improved: it=%d k=%d obj=%.6g\n', it, k, best_cost);
                end

                % Intensify
                k = kmin;
            else
                % No improvement: restore best & diversify
                x      = x_best;
                alphaA = alphaA_b;
                f1     = f1_b;
                f2     = f2_b;
                k      = k + kstep;
            end
        end
    end
end

% ======================================================================
% SUBFUNCTIONS
% ======================================================================

function val = get_opt(s, field, default)
    if isfield(s, field) && ~isempty(s.(field))
        val = s.(field);
    else
        val = default;
    end
end

function x = project_to_binary_k(x0, k_card)
    % Build a binary vector with exactly k_card ones.
    n = numel(x0);

    % Use x0 as score: larger entries are more likely to become 1.
    [~, perm] = sort(x0, 'descend');

    x = zeros(n,1);
    if k_card > 0
        x(perm(1:k_card)) = 1;
    end
end

function obj = cur_obj(f1, f2)
    obj = f1 + f2;
end

function [f1p, f2p] = trial_swap(i, j, x, f1, f2, alphaA, A, a)
    if x(i) ~= 1 || x(j) ~= 0
        error('trial_swap requires x(i)=1 and x(j)=0.');
    end

    f1p = f1 ...
        + 2 * (-alphaA(i) + alphaA(j)) ...
        + A(i,i) + A(j,j) - 2 * A(i,j);

    f2p = f2 - a(i) + a(j);
end

function [x, f1, f2, alphaA] = apply_swap(i, j, x, f1, f2, alphaA, A, a)
    if x(i) ~= 1 || x(j) ~= 0
        error('apply_swap requires x(i)=1 and x(j)=0.');
    end

    f1 = f1 ...
        + 2 * (-alphaA(i) + alphaA(j)) ...
        + A(i,i) + A(j,j) - 2 * A(i,j);

    f2 = f2 - a(i) + a(j);

    alphaA = alphaA - A(:,i) + A(:,j);

    x(i) = 0;
    x(j) = 1;
end

% Best-improvement local search in swap neighborhood
function [x, f1, f2, alphaA] = local_search_swap(x, f1, f2, alphaA, A, a, max_iters, max_time_sec, t0)
    EPS = 1e-4;
    it_count = 0;
    improved = true;

    while improved && it_count < max_iters
        if toc(t0) > max_time_sec
            return;
        end

        it_count = it_count + 1;
        improved = false;

        best_val = cur_obj(f1, f2);
        best_i   = 0;
        best_j   = 0;

        ones_idx  = find(x ~= 0);
        zeros_idx = find(x == 0);

        if isempty(ones_idx) || isempty(zeros_idx)
            return;
        end

        for p = 1:numel(ones_idx)
            if toc(t0) > max_time_sec
                return;
            end

            i = ones_idx(p);
            for q = 1:numel(zeros_idx)
                j = zeros_idx(q);

                [f1p, f2p] = trial_swap(i, j, x, f1, f2, alphaA, A, a);
                val = f1p + f2p;

                if val + EPS < best_val
                    best_val = val;
                    best_i   = i;
                    best_j   = j;
                end
            end
        end

        if best_i > 0
            [x, f1, f2, alphaA] = apply_swap(best_i, best_j, x, f1, f2, alphaA, A, a);
            improved = true;
        end
    end
end

% Shake: perform k random feasible swaps
function [x, f1, f2, alphaA] = shake(x, f1, f2, alphaA, A, a, k)
    for t = 1:k
        ones_idx  = find(x ~= 0);
        zeros_idx = find(x == 0);

        if isempty(ones_idx) || isempty(zeros_idx)
            return;
        end

        i = ones_idx(randi(numel(ones_idx)));
        j = zeros_idx(randi(numel(zeros_idx)));

        [x, f1, f2, alphaA] = apply_swap(i, j, x, f1, f2, alphaA, A, a);
    end
end