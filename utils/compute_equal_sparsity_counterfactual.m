function T=compute_equal_sparsity_counterfactual(Q,n_train,Ag,E,rel_tol,n_random)
% Original Figure 16 support selection, refit, and held-out evaluation.
Qtr=Q(1:n_train,:); Qfu=Q(n_train+1:end,:); r=size(Qtr,2);
[I,J]=ndgrid(1:r,1:r); off=I~=J; idxoff=find(off);
voff=abs(Ag(off)); thr=rel_tol*max(max(voff),eps); gcOff=(abs(Ag)>thr)&off; K=nnz(gcOff);
if K<1, error('GC effective support is empty at the requested threshold.'); end
K=min(K,numel(idxoff));

% All methods keep self-memory (diagonal) to match the GC penalty design.
base=eye(r)>0;
Sgc=base|gcOff;
sep=abs(I-J); [~,ord]=sort(sep(idxoff),'ascend'); Snear=base; Snear(idxoff(ord(1:K)))=true;
E=max(E(:),realmin); gap=abs(log10(E(I)./E(J))); [~,ord]=sort(gap(idxoff),'ascend'); Senergy=base; Senergy(idxoff(ord(1:K)))=true;

A_gc_refit=refit_support(Qtr,Sgc); A_near=refit_support(Qtr,Snear); A_energy=refit_support(Qtr,Senergy);
e_gc=heldout_error(A_gc_refit,Qtr,Qfu); e_near=heldout_error(A_near,Qtr,Qfu); e_energy=heldout_error(A_energy,Qtr,Qfu);

er=zeros(n_random,1); rng(42);
for b=1:n_random
  pick=idxoff(randperm(numel(idxoff),K)); Sr=base; Sr(pick)=true;
  Ar=refit_support(Qtr,Sr); er(b)=heldout_error(Ar,Qtr,Qfu);
end

% Also report actual regularized GC-DMD, but do not mix it into the equal-
% sparsity bars because its coefficients are shrinkage estimates, not LS refits.
e_actual=heldout_error(Ag,Qtr,Qfu);
active_ratio=K/numel(idxoff);

names={'GC-selected (refit)';'Near-diagonal';'Energy-local';'Random'};
means=[e_gc;e_near;e_energy;mean(er)]; stds=[0;0;0;std(er)];
T=table(names,means,stds,repmat(K,4,1),repmat(active_ratio,4,1),...
    'VariableNames',{'support_strategy','heldout_RMSE_mean','heldout_RMSE_std','active_offdiag_count','active_offdiag_ratio'});
T.actual_regularized_GC_RMSE=repmat(e_actual,4,1);
end

function A=refit_support(Qtr,S)
X=Qtr(1:end-1,:)'; Y=Qtr(2:end,:)'; r=size(X,1); A=zeros(r);
for i=1:r
  js=find(S(i,:));
  if isempty(js), continue; end
  Xi=X(js,:);
  % Minimum-norm constrained least-squares fit on the selected parents.
  A(i,js)=Y(i,:)*pinv(Xi);
end
end

function e=heldout_error(A,Qtr,Qfu)
% Paper definition: RMSE_test(A) = ||Q2_test - A*Q1_test||_F / sqrt(r*T_test).
% Use the final identification state to predict the first held-out state.
Qseq = [Qtr(end,:); Qfu];
Q1 = Qseq(1:end-1,:)';
Q2 = Qseq(2:end,:)';
residual = Q2 - A*Q1;
[r,T_test] = size(residual);
e = norm(residual,'fro') / sqrt(r*T_test);
end
