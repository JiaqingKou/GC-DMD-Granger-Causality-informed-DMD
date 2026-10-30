function [Sig, pcrit] = bh_fdr_mask(Pval, q, test_mask)
%BH_FDR_MASK Benjamini-Hochberg FDR correction for a matrix of p-values.
%
%   [Sig,pcrit] = bh_fdr_mask(Pval,q,test_mask)
%
% Pval      : matrix of p-values
% q         : target false-discovery rate (e.g. 0.05)
% test_mask : logical matrix selecting the family of hypotheses to correct
%
% Sig has the same size as Pval. Entries outside test_mask are false.

if nargin < 3 || isempty(test_mask)
    test_mask = true(size(Pval));
end
if q <= 0 || q >= 1
    error('q must lie in (0,1).');
end

pv = Pval(test_mask);
pv = pv(:);
valid = isfinite(pv);
pv_valid = pv(valid);

Sig = false(size(Pval));
pcrit = NaN;
if isempty(pv_valid)
    return;
end

[p_sorted,~] = sort(pv_valid,'ascend');
m = numel(p_sorted);
thresholds = q * (1:m)' / m;
k = find(p_sorted <= thresholds, 1, 'last');
if isempty(k)
    return;
end

pcrit = p_sorted(k);
Sig(test_mask) = Pval(test_mask) <= pcrit;
end
