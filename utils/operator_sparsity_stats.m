function S = operator_sparsity_stats(A,rel_tol)
%OPERATOR_SPARSITY_STATS Off-diagonal sparsity diagnostics.
r=size(A,1); mask=~eye(r); v=abs(A(mask));
if isempty(v), scale=0; else, scale=max(v); end
thr=rel_tol*max(scale,eps);
S.exact_zeros=nnz(v==0);
S.exact_zero_fraction=S.exact_zeros/max(numel(v),1);
S.active_offdiag=nnz(v>thr);
S.density=S.active_offdiag/max(numel(v),1);
S.effective_zero_fraction=1-S.density;
S.threshold=thr;
end
