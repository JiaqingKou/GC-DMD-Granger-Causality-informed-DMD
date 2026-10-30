function varargout = optimize_operator(mode,varargin)
%OPTIMIZE_OPERATOR Unified FISTA solver and validation-based lambda search.
%   [A,info] = optimize_operator('solve',data,lambda,W,opts)
%   R = optimize_operator('gc_path',train,val,W,grid,A_dmd,A_true,tol,opts,Sig,settings)
%   [lambda,A,rmse] = optimize_operator('validation',train,val,W,grid,opts)
% Modes share the SAME original numerical routines; GC weights are refreshed
% only for the gc_path final refit (and optional validation gc_settings).
switch lower(mode)
    case 'solve'
        [varargout{1:nargout}] = solve_one(varargin{:});
    case 'gc_path'
        [varargout{1:nargout}] = scan_gc(varargin{:});
    case 'validation'
        [varargout{1:nargout}] = scan_validation(varargin{:});
    otherwise
        error('Unknown optimization mode: %s',mode);
end
end

function [A, info] = solve_one(data, lambda, W, opts)
%FIT_WEIGHTED_L1_DMD Weighted-L1 regularized DMD using FISTA.
%
% Solve
%
%   min_A  (1/(2T))*||X2 - A X1||_F^2
%          + lambda * sum_ij W_ij |A_ij|
%
% Inputs:
%   data   : N x r
%   lambda : nonnegative scalar
%   W      : r x r nonnegative penalty weights
%   opts.max_iter
%   opts.tol
%   opts.verbose
%
% W_ii can be set to zero to leave self-memory unpenalized.

if nargin < 4
    opts = struct;
end
if ~isfield(opts,'max_iter'), opts.max_iter = 5000; end
if ~isfield(opts,'tol'), opts.tol = 1e-8; end
if ~isfield(opts,'verbose'), opts.verbose = false; end

X = data(1:end-1,:)';   % r x T
Y = data(2:end,:)';     % r x T
T = size(X,2);
r = size(X,1);

if ~isequal(size(W),[r r])
    error('W must be r x r.');
end
if any(W(:) < 0)
    error('W must be nonnegative.');
end
if lambda < 0
    error('lambda must be nonnegative.');
end

% Unregularized initialization
A0 = Y * pinv(X);

if lambda == 0
    A = A0;
    info.iter = 0;
    info.converged = true;
    info.objective = objective(A,Y,X,T,lambda,W);
    return;
end

% Lipschitz constant of gradient:
% grad f(A) = (A X - Y) X' / T
G = (X*X') / T;
Lips = norm(G,2);
if Lips <= eps
    error('Degenerate X matrix.');
end
step = 1 / Lips;

A = A0;
Z = A;
t = 1;

obj_prev = objective(A,Y,X,T,lambda,W);
converged = false;

for it = 1:opts.max_iter

    grad = ((Z*X - Y) * X') / T;

    V = Z - step*grad;
    thresh = step * lambda .* W;

    Anew = soft_threshold(V, thresh);

    tnew = 0.5*(1 + sqrt(1 + 4*t^2));
    Znew = Anew + ((t - 1)/tnew)*(Anew - A);

    relchg = norm(Anew-A,'fro') / max(norm(A,'fro'),1e-12);

    A = Anew;
    Z = Znew;
    t = tnew;

    if mod(it,25)==0 || it==1
        obj = objective(A,Y,X,T,lambda,W);

        if opts.verbose
            fprintf('  it=%5d  obj=%.6e  relchg=%.3e\n', it, obj, relchg);
        end

        % Monotonicity safeguard: if accelerated step overshoots strongly,
        % restart FISTA momentum.
        if obj > obj_prev*(1 + 1e-10)
            Z = A;
            t = 1;
        end
        obj_prev = obj;
    end

    if relchg < opts.tol
        converged = true;
        break;
    end
end

info.iter = it;
info.converged = converged;
info.objective = objective(A,Y,X,T,lambda,W);
info.step = step;
info.Lipschitz = Lips;

end


function Z = soft_threshold(X, tau)
Z = sign(X) .* max(abs(X)-tau,0);
end


function obj = objective(A,Y,X,T,lambda,W)
R = Y - A*X;
obj = 0.5/T * sum(R(:).^2) + lambda * sum((W.*abs(A)),'all');
end


function R = scan_gc(dataTrain,dataVal,W,lambda_grid,A_dmd,A_true,tol,opts,GC_support,gc_settings)
%ANALYZE_GC_LAMBDA_PATH Track GC-DMD operator A(lambda) along lambda.
%
% Optional GC_support is the binary significant GC graph (j -> i stored at i,j).
% If supplied, the routine also measures how closely supp(A_lambda) matches
% the significant GC graph using Precision / Recall / F1 / Jaccard.
if nargin < 7 || isempty(tol), tol = 1e-3; end
if nargin < 9, GC_support = []; end

if nargin<10 || isempty(gc_settings), error('GC refit settings are required.'); end
nL = numel(lambda_grid); r = size(A_dmd,1); off = ~eye(r);
R.lambda = lambda_grid(:);
R.A = cell(nL,1);
R.val_rmse = zeros(nL,1);
R.rel_from_dmd = zeros(nL,1);
R.rel_to_true = zeros(nL,1);
R.offdiag_density = zeros(nL,1);
R.spectral_radius = zeros(nL,1);
R.l1_offdiag = zeros(nL,1);

if ~isempty(GC_support)
    G = logical(GC_support) & off;
    R.gc_precision = zeros(nL,1);
    R.gc_recall    = zeros(nL,1);
    R.gc_f1        = zeros(nL,1);
    R.gc_jaccard   = zeros(nL,1);
    R.gc_edge_count = nnz(G);
else
    G = [];
end

for k = 1:nL
    [A,~] = solve_one(dataTrain,lambda_grid(k),W,opts);
    R.A{k} = A;
    R.val_rmse(k) = one_step_rmse(A,dataVal);
    R.rel_from_dmd(k) = norm(A-A_dmd,'fro')/max(norm(A_dmd,'fro'),eps);
    R.rel_to_true(k) = norm(A-A_true,'fro')/max(norm(A_true,'fro'),eps);
    Asupp = (abs(A)>tol) & off;
    R.offdiag_density(k) = nnz(Asupp)/nnz(off);
    R.spectral_radius(k) = max(abs(eig(A)));
    R.l1_offdiag(k) = sum(abs(A(off)));

    if ~isempty(G)
        TP = nnz(Asupp & G);
        FP = nnz(Asupp & ~G & off);
        FN = nnz(~Asupp & G);
        P = TP/max(TP+FP,1);
        Rc = TP/max(TP+FN,1);
        F1 = 2*P*Rc/max(P+Rc,eps);
        J = TP/max(TP+FP+FN,1);
        R.gc_precision(k)=P;
        R.gc_recall(k)=Rc;
        R.gc_f1(k)=F1;
        R.gc_jaccard(k)=J;
    end
end

[~,R.best_idx] = min(R.val_rmse);
R.best_lambda = R.lambda(R.best_idx);
R.best_A_train = R.A{R.best_idx};
% Validation chooses lambda; final identification uses train + validation.
[Wfinal,R.GC_raw_final,R.Sig_final,R.C_gc_final]=build_gc_weights([dataTrain;dataVal],gc_settings);
R.W_final=Wfinal;
[R.best_A,R.refit_info] = solve_one([dataTrain;dataVal],R.best_lambda,Wfinal,opts);
R.refit_after_validation = true;
R.refit_snapshot_count = size(dataTrain,1)+size(dataVal,1);
R.support_tol = tol;

if ~isempty(G)
    maxF = max(R.gc_f1);
    cand = find(abs(R.gc_f1-maxF) <= 1e-12);
    % Tie-break: among equally good support matches, prefer the one nearest
    % the validation-selected lambda on the log scale.
    if numel(cand)>1
        lp = log10(max(R.lambda(cand),realmin));
        lb = log10(max(R.best_lambda,realmin));
        [~,ii] = min(abs(lp-lb));
        cand = cand(ii);
    end
    R.gc_best_idx = cand(1);
    R.gc_best_lambda = R.lambda(R.gc_best_idx);
    R.gc_best_f1 = R.gc_f1(R.gc_best_idx);
    R.gc_support = G;
    R.gc_best_A_train=R.A{R.gc_best_idx};
    R.gc_best_A=solve_one([dataTrain;dataVal],R.gc_best_lambda,Wfinal,opts);
end
end


function [best_lambda, best_A, val_rmse] = scan_validation( ...
    dataTrain, dataVal, W, lambda_grid, fista_opts, gc_settings)
%SELECT_LAMBDA_VALIDATION Select lambda using chronological validation.
%
% Selection metric: one-step RMSE on the validation sequence.

nL = numel(lambda_grid);
val_rmse = zeros(nL,1);
Acell = cell(nL,1);

for k = 1:nL
    lam = lambda_grid(k);

    [A,~] = solve_one(dataTrain, lam, W, fista_opts);

    val_rmse(k) = one_step_rmse(A, dataVal);
    Acell{k} = A;
end

[~,idx] = min(val_rmse);

best_lambda = lambda_grid(idx);
% Optional GC refresh; Oracle weights remain prescribed.
if nargin>=6 && ~isempty(gc_settings)
 W=build_gc_weights([dataTrain;dataVal],gc_settings);
end
[best_A,~] = solve_one([dataTrain;dataVal],best_lambda,W,fista_opts);

end
