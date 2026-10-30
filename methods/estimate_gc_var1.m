function [GC_raw, Pval, Sig, C, info] = estimate_gc_var1(data, varargin)
%ESTIMATE_GC_VAR1 Classical conditional multivariate lag-1 Granger causality.
%
% data : N x r
% Convention: GC_raw(i,j) = j -> i
%
% For each target i:
%   full model    : x_i(t+1) ~ 1 + all x(t)
%   reduced model : x_i(t+1) ~ 1 + all x(t) except x_j(t)
%   The intercept is included in BOTH models and is never removed.
%   GC uses one-step OLS residuals; no recursive prediction is used.
%
% GC strength:
%   F_GC(j->i) = log( var(e_red) / var(e_full) )
%
% Significance:
%   standard nested linear-model F test, followed by optional BH-FDR
%   correction across all off-diagonal directed edges.
%
% Options:
%   'alpha'            : raw per-test level, default 0.05
%   'multiple_testing' : 'bh-fdr' (default) or 'none'
%   'fdr_q'            : BH target FDR, default = alpha
%   'norm_mode'        : 'rowsum' (default) or 'rowmax'
%   'intercept'        : true (default); include a constant in both regressions

p = inputParser;
addRequired(p,'data');
addParameter(p,'alpha',0.05);
addParameter(p,'multiple_testing','bh-fdr');
addParameter(p,'fdr_q',[]);
addParameter(p,'norm_mode','rowsum');
addParameter(p,'intercept',true);
parse(p,data,varargin{:});

alpha = p.Results.alpha;
mt_mode = lower(string(p.Results.multiple_testing));
fdr_q = p.Results.fdr_q;
if isempty(fdr_q), fdr_q = alpha; end
norm_mode = lower(p.Results.norm_mode);
use_intercept = logical(p.Results.intercept);

[N,r] = size(data);
if N <= r + 3
    error('Too few samples for full VAR(1) with intercept: N must exceed r+3.');
end

X = data(1:end-1,:);
Y = data(2:end,:);
T = size(X,1);

% Design matrix: intercept + r lag-one predictors.
% Do not penalize/remove the intercept when testing j -> i.
if use_intercept
    X_full = [ones(T,1), X];
else
    X_full = X;
end
B_full = X_full \ Y;
Yhat_full = X_full * B_full;
E_full = Y - Yhat_full;
RSS_full = sum(E_full.^2,1);
var_full = RSS_full / T;

GC_raw = zeros(r);
Pval   = ones(r);
Sig_raw = false(r);

K_full = r + double(use_intercept); % includes constant coefficient
df2 = T - K_full;

for i = 1:r
    y = Y(:,i);
    for j = 1:r
        if i == j, continue; end

        keep = true(1,r); keep(j) = false;
        Xred = X(:,keep);
        if use_intercept
            Xred = [ones(T,1), Xred]; % same intercept in reduced model
        end
        b_red = Xred \ y;
        e_red = y - Xred*b_red;

        RSS_red = sum(e_red.^2);
        var_red = RSS_red / T;
        gc = log(max(var_red,eps) / max(var_full(i),eps));
        GC_raw(i,j) = max(gc,0);

        df1 = 1; % exactly one lagged predictor removed; intercept retained
        numerator = max(RSS_red - RSS_full(i), 0) / df1;
        denominator = max(RSS_full(i),eps) / max(df2,1);
        Fstat = numerator / denominator;

        z = (df1 * Fstat) / (df1 * Fstat + df2);
        z = min(max(z,0),1);
        Fcdf = betainc(z, df1/2, df2/2);
        pval = 1 - Fcdf;

        Pval(i,j) = pval;
        Sig_raw(i,j) = (pval < alpha);
    end
end

% Multiple-testing correction over the complete directed off-diagonal family.
offdiag = ~eye(r);
switch mt_mode
    case {"bh-fdr","bh","fdr"}
        [Sig, pcrit] = bh_fdr_mask(Pval, fdr_q, offdiag);
    case {"none","raw"}
        Sig = Sig_raw;
        pcrit = alpha;
    otherwise
        error('Unknown multiple_testing mode: %s', mt_mode);
end
Sig(1:r+1:end) = false;

GC_sig = GC_raw .* Sig;
GC_sig(1:r+1:end) = 0;

C = zeros(r);
switch norm_mode
    case 'rowsum'
        for i = 1:r
            s = sum(GC_sig(i,:));
            if s > 0, C(i,:) = GC_sig(i,:) / s; end
        end
    case 'rowmax'
        for i = 1:r
            m = max(GC_sig(i,:));
            if m > 0, C(i,:) = GC_sig(i,:) / m; end
        end
    otherwise
        error('Unknown norm_mode: %s', norm_mode);
end
C(1:r+1:end) = 0;

info.B_full = B_full;
info.intercept = use_intercept;
if use_intercept
    info.b_full = B_full(1,:); % estimated intercept for each target
    info.A_full = B_full(2:end,:).'; % rows: effects, columns: causes
else
    info.b_full = zeros(1,r);
    info.A_full = B_full.';
end
info.E_full = E_full;
info.var_full = var_full;
info.GC_sig = GC_sig;
info.Sig_raw = Sig_raw;
info.alpha = alpha;
info.multiple_testing = char(mt_mode);
info.fdr_q = fdr_q;
info.pcrit = pcrit;
info.norm_mode = norm_mode;
info.T = T;
info.df2 = df2;
end
