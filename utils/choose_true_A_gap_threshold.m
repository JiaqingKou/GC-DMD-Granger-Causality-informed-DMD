function tau=choose_true_A_gap_threshold(A)
% Ground-truth-only cutoff: midpoint of the gap between exact zeros and
% smallest nonzero off-diagonal magnitude. No GC information is used.
r=size(A,1); vals=abs(A(~eye(r))); nz=vals(vals>0);
if isempty(nz), error('True A has no nonzero off-diagonal entries.'); end
if ~any(vals==0), error('True A has no exact-zero off-diagonal entries; no zero gap exists.'); end
tau=min(nz)/2;
fprintf('True-A exact-zero gap: max zero=0, min nonzero=%.8g, midpoint threshold=%.8g\n',min(nz),tau);
end
