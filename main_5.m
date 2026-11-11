clear; clc; close all; rng(42);
%% ---------------- User settings ----------------
root=fileparts(mfilename('fullpath'));
addpath(fullfile(root,'methods'),fullfile(root,'utils'),fullfile(root,'figures'));
data_file=fullfile(root,'data','Cavity30k.mat');
output_dir=prepare_results_dir(root,5);
error_curve_ylim=[0 6]; validation_ylim=[0 .012]; % paper display limits
pod_energy_target=0.99999;
gc_color_max=[]; % []: automatic off-diagonal positive 99th percentile; or explicit maximum
operator_color_max=0.05; % shared operator color limits +/-0.05
operator_zero_tol_rel=1e-3;
n_model=800; n_fit=640; n_val=160; pred_horizon=400;
% Spatial-field snapshots are set here (absolute snapshot indices).
% Reconstruction: choose indices inside 1:n_model. Prediction: choose indices inside n_model+1:n_model+pred_horizon.
reconstruction_field_snapshots=800;
prediction_field_snapshots=[1000 1200];
% Three Appendix B field snapshots (B2b–B2d).
error_field_snapshots=[800 1000 1200];
n_lambda=33;
fista_opts.max_iter=5000; fista_opts.tol=1e-8; fista_opts.verbose=false;
fprintf('\n=========== GC-DMD Re=30000 cavity benchmark ===========\n');

%% Step 1. Load snapshots
D=load_cavity_snapshots(data_file,'VelocityField','Psi',n_model,pred_horizon);
Xraw=D.Xraw; N=D.N; time=D.time; dt=D.dt;
mu=mean(Xraw(:,1:n_model),2);
Xc=Xraw-mu; % Mean subtraction is performed once, here in main.

%% 2. SVD of mean-subtracted snapshots
[U,sig,V]=compute_svd(Xc(:,1:n_model));
cum_energy=cumsum(sig.^2)/sum(sig.^2);
r=find(cum_energy>=pod_energy_target,1);
[Ur,~]=qr(U(:,1:r),0);
Q=(Ur'*Xc)'; Qmodel=Q(1:n_model,:);
fprintf('POD rank r=%d, retained energy=%.8f\n',r,cum_energy(r));

%% 3. Standard DMD: least-squares identification on reduced coordinates
% The first 640 modeling snapshots form the fitting subset; the next 160
% are held out for validation. The independent prediction horizon is not
% used for operator fitting or regularization-parameter selection.
train=Qmodel(1:n_fit,:);
val=Qmodel(n_fit+1:end,:);
% Stage 1: fit DMD on the training subset only. This temporary operator
% provides a baseline one-step RMSE on the same validation snapshots used
% to select the GC-DMD regularization parameter; it is not the final DMD.
Xfit=train(1:end-1,:)'; Yfit=train(2:end,:)';
A_dmd_fit=Yfit*pinv(Xfit);
val_dmd=one_step_rmse(A_dmd_fit,val);
% Stage 2: after validation, refit standard DMD on all modeling snapshots
% (training + validation). This is the operator used for final comparisons.
X1=Qmodel(1:end-1,:)'; X2=Qmodel(2:end,:)';
A_dmd=X2*pinv(X1); % final fit on all 800 modeling snapshots

%% 4. Conditional Granger causality: full/reduced VAR(1) with intercept
% During parameter selection, infer conditional GC from training data only.
% Its BH-FDR-filtered, row-normalized structure sets the relative penalty
% weights; the validation subset is not used to estimate these weights.
gc_args={'alpha',0.05,'multiple_testing','bh-fdr','fdr_q',0.05, ...
         'norm_mode','rowsum','intercept',true};
[~,~,~,C_train]=estimate_gc_var1(train,gc_args{:});
W_train=build_gc_weights(C_train); % eta=1, diagonal unpenalized

%% 5. GC-DMD: choose lambda_pred, then identify the final operator
% Select lambda using validation one-step RMSE for operators fitted on
% training data; then reestimate GC on training + validation and refit
% GC-DMD with the SAME selected lambda (no test-data leakage).
lambda_grid=build_lambda_grid(train,n_lambda,'cavity30k');
[lambda_gc,~,val_gc]=optimize_operator('validation',train,val,W_train,lambda_grid,fista_opts);
[GC_raw,Pval,Sig,C_gc_final,gc_info]=estimate_gc_var1(Qmodel,gc_args{:}); %#ok<NASGU>
W_gc=build_gc_weights(C_gc_final);
A_gc=optimize_operator('solve',Qmodel,lambda_gc,W_gc,fista_opts);
val=Qmodel(n_fit+1:n_fit+n_val,:); % held-out validation coordinates

%% 6. Operator diagnostics
fprintf('DMD and GC-DMD operators identified; lambda_pred=%.6g\n',lambda_gc);

%% 7. Recursive prediction and vorticity fields for Appendix B Figure B2
figureB2_indices=unique([reconstruction_field_snapshots prediction_field_snapshots error_field_snapshots],'stable');
% Starting from the same observed reduced state, repeatedly apply each
% final operator without injecting subsequent true states (free-run forecast).
% This differs from one-step validation, which uses observed previous states.
nt=min(size(Xraw,2),n_model+pred_horizon);
q0=Q(1,:)';
qd=zeros(nt,size(Q,2)); qg=qd;
qd(1,:)=q0'; qg(1,:)=q0';
for k=2:nt
    qd(k,:)=(A_dmd*qd(k-1,:)')';
    qg(k,:)=(A_gc*qg(k-1,:)')';
end
% Reconstruct full fields and evaluate normalized snapshot errors.
Xd=mu+Ur*qd'; Xg=mu+Ur*qg'; Xtruth=Xraw(:,1:nt);
den=sqrt(sum((Xtruth-mu).^2,1)); den=max(den,eps);
ed=sqrt(sum((Xd-Xtruth).^2,1))./den;
eg=sqrt(sum((Xg-Xtruth).^2,1))./den;
figureB2_data=struct('nt',nt,'ed',ed,'eg',eg,'Xtruth',Xtruth,'Xd',Xd,'Xg',Xg);
% The helper below only prepares spatial/vorticity fields for plotting;
% time marching, reconstruction, and error evaluation remain visible here.
figureB2_data=compute_cavity_figure_data(figureB2_data,N,figureB2_indices,'cavity30k');

%% 7. Appendix B Figure B1 (original 30k four-panel layout)
figB1=plot_figB1(C_gc_final,lambda_grid,val_gc,val_dmd,lambda_gc, ...
    A_dmd,A_gc,operator_color_max,gc_color_max,'',validation_ylim);
print(figB1,fullfile(output_dir,'B1.eps'),'-depsc2','-painters');

%% 8. Appendix B Figure B2 (original 30k curve and three field windows)
[figB2a,figB2b,figB2c,figB2d]=plot_figB2(figureB2_data,N, ...
    n_model,pred_horizon,dt,reconstruction_field_snapshots, ...
    prediction_field_snapshots,[],error_field_snapshots,time,false,error_curve_ylim);
print(figB2a,fullfile(output_dir,'B2a.png'),'-dpng','-r300');
print(figB2b,fullfile(output_dir,'B2b.png'),'-dpng','-r300');
print(figB2c,fullfile(output_dir,'B2c.png'),'-dpng','-r300');
print(figB2d,fullfile(output_dir,'B2d.png'),'-dpng','-r300');
fprintf('\n30k Appendix B figures exported to %s\n',output_dir);
