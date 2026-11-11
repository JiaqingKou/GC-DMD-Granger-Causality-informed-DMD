function R=identify_operators(Qmodel,n_fit,lambda_grid,opts)
% Common cavity/jet: train -> validation lambda -> full-model GC refit.
train=Qmodel(1:n_fit,:);
val=Qmodel(n_fit+1:end,:);
[~,~,~,~,Wtrain]=infer_gc_weights(train);
Xtr=train(1:end-1,:)'; Ytr=train(2:end,:)'; A_dmd_fit=Ytr*pinv(Xtr);
val_dmd=one_step_rmse(A_dmd_fit,val);
[lambda_gc,~,val_gc]=optimize_operator('validation',train,val,Wtrain,lambda_grid,opts);
[GC_raw,Pval,Sig,C_gc,Wgc,info]=infer_gc_weights(Qmodel);
X1=Qmodel(1:end-1,:)'; X2=Qmodel(2:end,:)'; A_dmd=X2*pinv(X1);
A_gc=optimize_operator('solve',Qmodel,lambda_gc,Wgc,opts);
R=struct('A_dmd',A_dmd,'A_gc',A_gc,'A_dmd_fit',A_dmd_fit, ...
    'lambda_gc',lambda_gc,'lambda_grid',lambda_grid,'val_gc',val_gc, ...
    'cv_gc',val_gc,'cv_dmd',val_dmd,'GC_raw',GC_raw,'Pval',Pval, ...
    'Sig',Sig,'Sig_final',Sig,'C_gc',C_gc,'C_gc_final',C_gc, ...
    'gc_info',info);
end
