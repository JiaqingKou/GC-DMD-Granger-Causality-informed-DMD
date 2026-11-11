clear; clc; close all; rng(42);
%% ---------------- User settings ----------------
root=fileparts(mfilename('fullpath'));
addpath(fullfile(root,'methods'),fullfile(root,'utils'),fullfile(root,'figures'));
data_file=fullfile(root,'data','jetLES.mat');
output_dir=prepare_results_dir(root,6);
%% 0. Parameters for the Jet benchmark (paper Figures 13--16)
cfg.n_train=640; cfg.n_future=160; cfg.pod_rank=300;
cfg.rank_list=[100 125 150 200 250 300];
cfg.alpha_gc=.05; cfg.gc_norm_mode='rowsum'; cfg.gc_eta=1;
cfg.val_frac_inside_train=.20; cfg.n_lambda=67; cfg.lambda_exp_range=[-10 1]; cfg.pareto_max_points=15;
cfg.fista.max_iter=5000; cfg.fista.tol=1e-8; cfg.fista.verbose=false;
cfg.support_tol=1e-4; cfg.dt=.2;
%% Step 1. Load pressure snapshots
fn=data_file;
if ~isfile(fn),error('Copy jetLES.mat into data/ next to methods/.');end
S=load(fn); if ~isfield(S,'p'),error('jetLES.mat must contain p.');end
nx=39; ny=175; ns=5000;
Xraw=(reshape(double(S.p),ns,nx*ny))';
nuse=cfg.n_train+cfg.n_future;
mu=mean(Xraw(:,1:cfg.n_train),2);
Xc=bsxfun(@minus,Xraw(:,1:nuse),mu); % center immediately after import
%% Step 2. SVD/POD basis and reduced coefficients
% Construct the basis from the 640 modeling snapshots only; project the
% subsequent held-out snapshots onto this fixed basis for evaluation.
[U,s,~]=compute_svd(Xc(:,1:cfg.n_train));
basis=struct('U',U(:,1:cfg.pod_rank),'singular_values',s, ...
    'energy',s.^2/sum(s.^2),'cum_energy',cumsum(s.^2)/sum(s.^2), ...
    'rank',cfg.pod_rank);
Ur=basis.U;
Q=(Ur'*Xc)';

%% Step 3. Standard DMD on the reduced training coordinates
% Two DMD fits serve different purposes: the first uses only the fitting
% subset to obtain a comparable validation RMSE; the second uses all 640
% modeling snapshots (fit + validation) for the final DMD operator.
% Reserve the final 20%% of the modeling interval for validation.
% The rank-dependent lower bound ensures enough fitting snapshot pairs.
n_fit=floor((1-cfg.val_frac_inside_train)*cfg.n_train);
n_fit=max(n_fit,cfg.pod_rank+3); n_fit=min(n_fit,cfg.n_train-2);
Qmodel=Q(1:cfg.n_train,:);
train=Qmodel(1:n_fit,:); val=Qmodel(n_fit+1:end,:);
Xfit=train(1:end-1,:)'; Yfit=train(2:end,:)';
% Temporary training-only DMD, evaluated on unseen validation snapshots.
A_dmd_fit=Yfit*pinv(Xfit);
val_dmd=one_step_rmse(A_dmd_fit,val);
% Final DMD refit on the entire modeling interval, not the future test set.
X1=Qmodel(1:end-1,:)'; X2=Qmodel(2:end,:)';
A_dmd=X2*pinv(X1);

%% Step 4. Conditional GC with intercept and BH-FDR
% Training-only conditional GC supplies the off-diagonal penalty weights
% for lambda selection; no validation or future snapshots enter this step.
% The final GC structure will be reestimated after lambda selection.
gc_args={'alpha',0.05,'multiple_testing','bh-fdr','fdr_q',0.05, ...
         'norm_mode','rowsum','intercept',true};
[~,~,~,C_train]=estimate_gc_var1(train,gc_args{:});
W_train=build_gc_weights(C_train);

%% Step 5. GC-DMD: validation-selected regularization and final refit
% Choose lambda from validation one-step RMSE, then recompute GC using
% all modeling snapshots and solve the final GC-DMD operator with that lambda.
lambda_grid=build_lambda_grid(train,cfg.n_lambda,'jet');
[lambda_gc,~,val_gc]=optimize_operator('validation',train,val,W_train,lambda_grid,cfg.fista);
[GC_raw,Pval,Sig,C_gc,gc_info]=estimate_gc_var1(Qmodel,gc_args{:});
W_gc=build_gc_weights(C_gc);
A_gc=optimize_operator('solve',Qmodel,lambda_gc,W_gc,cfg.fista);
fprintf('Jet DMD/GC-DMD identified; lambda_pred=%.6g\n',lambda_gc);
% The following numerical studies feed the manuscript figures; their
% plotting routines receive the computed results rather than fitting models.
%% 6. Rank-dependent operator sparsity (Figure 14)
Rrank=rank_sparsity_study(Xc,cfg);

%% 7. Validation and test sparsity--generalization (Figure 15)
P=build_gc_pareto(Q,cfg.n_train,cfg);

%% 8. Equal-sparsity counterfactual comparison (Figure 16)
cf_rank=size(A_gc,1);
cf_energy=basis.singular_values(1:cf_rank).^2;
cf_energy=cf_energy(:)/sum(cf_energy);
counterfactual=compute_equal_sparsity_counterfactual(Q,cfg.n_train,A_gc,cf_energy,cfg.support_tol,30);

%% 9. Paper figures: draw and export in manuscript order
fig13=plot_fig13(Xc,basis,nx,ny,[1 320 640 800]);
print(fig13,fullfile(output_dir,'fig13.png'),'-dpng','-r300');
fig14=plot_fig14(C_gc,A_dmd,A_gc,Rrank);
print(fig14,fullfile(output_dir,'fig14.eps'),'-depsc2','-painters');
% Figure 15 expects only four regularization-selection fields.
fig15_data=struct('lambda_grid',lambda_grid,'cv_gc',val_gc, ...
    'cv_dmd',val_dmd,'lambda_gc',lambda_gc);
[fig15a,fig15b]=plot_fig15(fig15_data,P);
print(fig15a,fullfile(output_dir,'fig15a.eps'),'-depsc2','-painters');
print(fig15b,fullfile(output_dir,'fig15b.eps'),'-depsc2','-painters');
fig16=plot_fig16(counterfactual);
print(fig16,fullfile(output_dir,'fig16.eps'),'-depsc2','-painters');
fprintf('Jet paper Figures 13--16 exported to %s\n',output_dir);
