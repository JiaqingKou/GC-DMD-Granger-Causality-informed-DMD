function R=rank_sparsity_study(X,cfg)
rlist=cfg.rank_list(:); nr=numel(rlist); act=nan(nr,2); err=nan(nr,2); energy=nan(nr,1); diagfrac=nan(nr,2);
band05=nan(nr,1); band10=nan(nr,1); band15=nan(nr,1);
Xtr=X(:,1:cfg.n_train); Xp=X(:,cfg.n_train+1:end);
[U,s,~]=compute_svd(Xtr);
cumulative=cumsum(s.^2)/sum(s.^2);
for k=1:nr
 r=rlist(k); fprintf('Rank sparsity study: r=%d\n',r); energy(k)=cumulative(r); Qtr=(U(:,1:r)'*Xtr)'; Qp=(U(:,1:r)'*Xp)'; nf=floor((1-cfg.val_frac_inside_train)*cfg.n_train); nf=max(nf,r+3); nf=min(nf,cfg.n_train-2); grid=build_lambda_grid(Qtr(1:nf,:),cfg.n_lambda,'jet'); Rk=identify_operators(Qtr,nf,grid,cfg.fista); A={Rk.A_dmd,Rk.A_gc};
 for m=1:2
   act(k,m)=nnz(abs(A{m})>cfg.support_tol)/r^2;
   diagfrac(k,m)=sum(abs(diag(A{m})))/sum(abs(A{m}),'all');
 end
 G=abs(Rk.A_gc); den=sum(G,'all'); [ii,jj]=ndgrid(1:r,1:r); sep=abs(ii-jj);
 band05(k)=sum(G(sep<=ceil(.05*r)),'all')/den;
 band10(k)=sum(G(sep<=ceil(.10*r)),'all')/den;
 band15(k)=sum(G(sep<=ceil(.15*r)),'all')/den;
 Qseq=[Qtr(end,:);Qp]; for m=1:2, E=(A{m}*Qseq(1:end-1,:)')'-Qseq(2:end,:); err(k,m)=mean(vecnorm(E,2,2)./max(vecnorm(Qseq(2:end,:),2,2),1e-12)); end
end
R=table(rlist,energy,act(:,1),act(:,2),diagfrac(:,1),diagfrac(:,2),band05,band10,band15,err(:,1),err(:,2), ...
 'VariableNames',{'rank','energy','active_DMD','active_GC','diagfrac_DMD','diagfrac_GC','band05_GC','band10_GC','band15_GC','error_DMD','error_GC'});
end
