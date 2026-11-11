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
%   standard nested linear-model F test, followed by BH-FDR
%   correction across all off-diagonal directed edges.
%
% Manuscript-fixed choices: intercept in both regressions, BH-FDR q=0.05,
% row-sum normalization, and exclusion of self-links.

% varargin is accepted for compatibility with existing main scripts;
% all manuscript cases use the same fixed BH-FDR/rowsum/intercept method.
alpha = 0.05;
[N,r] = size(data);
if N <= r + 3
    error('Too few samples for full VAR(1) with intercept: N must exceed r+3.');
end

X = data(1:end-1,:);
Y = data(2:end,:);
T = size(X,1);

% Design matrix: intercept + r lag-one predictors.
% Do not penalize/remove the intercept when testing j -> i.
X_full = [ones(T,1), X];
B_full = X_full \ Y;
Yhat_full = X_full * B_full;
E_full = Y - Yhat_full;
RSS_full = sum(E_full.^2,1);
var_full = RSS_full / T;

GC_raw = zeros(r);
Pval   = ones(r);
Sig_raw = false(r);

K_full = r + 1; % includes constant coefficient
df2 = T - K_full;

for i = 1:r
    y = Y(:,i);
    for j = 1:r
        if i == j, continue; end

        keep = true(1,r); keep(j) = false;
        Xred = X(:,keep);
        Xred = [ones(T,1), Xred]; % same intercept in reduced model
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

% Paper method: BH-FDR (q=0.05) on all directed off-diagonal tests.
offdiag = ~eye(r);
[Sig, pcrit] = bh_fdr_mask(Pval, 0.05, offdiag);
Sig(1:r+1:end) = false;
GC_sig = GC_raw .* Sig;
GC_sig(1:r+1:end) = 0;

% Paper method: row-sum normalization; zero rows remain zero.
row_sum = sum(GC_sig,2);
row_sum(row_sum==0) = 1;
C = bsxfun(@rdivide,GC_sig,row_sum);
C(1:r+1:end) = 0;

info.B_full = B_full;
info.intercept = true;
info.b_full = B_full(1,:);
info.A_full = B_full(2:end,:).';
info.E_full = E_full;
info.var_full = var_full;
info.GC_sig = GC_sig;
info.Sig_raw = Sig_raw;
info.alpha = alpha;
info.multiple_testing = 'bh-fdr';
info.fdr_q = 0.05;
info.pcrit = pcrit;
info.norm_mode = 'rowsum';
info.T = T;
info.df2 = df2;
end
