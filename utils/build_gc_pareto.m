function P=build_gc_pareto(Qall,n_train,cfg)
Qtr=Qall(1:n_train,:); Qp=Qall(n_train+1:end,:); r=size(Qtr,2); nTune=floor((1-cfg.val_frac_inside_train)*n_train); nTune=max(nTune,r+3); nTune=min(nTune,n_train-2); Qt=Qtr(1:nTune,:); Qv=Qtr(nTune+1:end,:);

[~,~,~,~,Wt]=infer_gc_weights(Qt); lg=build_lambda_grid(Qt,cfg.n_lambda,'jet'); [lam,~,~]=optimize_operator('validation',Qt,Qv,Wt,lg,cfg.fista);
[~,~,~,~,W]=infer_gc_weights(Qtr);
ids=unique(round(linspace(1,numel(lg),min(cfg.pareto_max_points,numel(lg))))); [~,ii]=min(abs(lg-lam)); ids=unique([ids(:);ii]); lams=lg(ids); n=numel(lams); active=zeros(n,1); pred=zeros(n,1);
for k=1:n, A=optimize_operator('solve',Qtr,lams(k),W,cfg.fista); active(k)=nnz(abs(A)>cfg.support_tol)/r^2; pred(k)=prmse(A,Qtr,Qp); end
X1=Qtr(1:end-1,:)'; X2=Qtr(2:end,:)'; Ad=X2*pinv(X1); P=struct('lambda',lams,'active_gc',active,'pred_gc',pred,'active_dmd',nnz(abs(Ad)>cfg.support_tol)/r^2,'pred_dmd',prmse(Ad,Qtr,Qp),'selected_lambda',lam,'selected_idx',find(ids==ii,1));
end
function e=prmse(A,Qtr,Qp), Q=[Qtr(end,:);Qp]; d=(A*Q(1:end-1,:)')'-Q(2:end,:); e=sqrt(mean(d(:).^2)); end
