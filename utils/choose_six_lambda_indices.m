function idx=choose_six_lambda_indices(R)
lam=R.lambda(:); n=numel(lam); ib=R.best_idx; ig=R.gc_best_idx;
positive=find(lam>0); before=find(lam>0 & lam<lam(ib));
if isempty(before), iw=positive(1); else
 target=sqrt(lam(positive(1))*max(lam(ib),lam(positive(1))));
 [~,k]=min(abs(log(lam(before))-log(target))); iw=before(k);
end
ref=max(lam(ib),lam(ig)); after=find(lam>ref & lam<max(lam));
if isempty(after), is=max(1,n-1); else
 target=sqrt(ref*max(lam)); [~,k]=min(abs(log(lam(after))-log(target))); is=after(k);
end
idx=[1 iw ib ig is n];
% Preserve six distinct stages if the optimum indices coincide.
for m=2:5
 if any(idx(m)==idx(1:m-1))
  available=setdiff(2:n-1,idx([1:m-1 m+1:end]));
  if ~isempty(available), [~,q]=min(abs(available-idx(m))); idx(m)=available(q); end
 end
end
end
