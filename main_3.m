clear; clc; close all; rng(42);
%% GC-DMD benchmark: Re=13000 lid-driven cavity flow
% Paper-consistent revision: DMD vs GC-DMD only; all figures follow project-1 style.

%% ---------------- User settings ----------------
root=fileparts(mfilename('fullpath'));
addpath(fullfile(root,'methods'),fullfile(root,'utils'),fullfile(root,'figures'));
data_file=fullfile(root,'data','Cavity13k.mat');
output_dir=prepare_results_dir(root,3); 
pod_energy_target=0.99999;
operator_color_max=0.05; % shared operator color limits +/-0.05
operator_zero_tol_rel=1e-3; modal_top_n=12;
pod_basis_modes=[1 3 5 7 9 11]; % one member from each leading POD pair
n_model=800; n_fit=640; n_val=160; pred_horizon=2200;
% Spatial-field snapshots are set here (absolute snapshot indices).
% Reconstruction: choose indices inside 1:n_model. Prediction: choose indices inside n_model+1:n_model+pred_horizon.
reconstruction_field_snapshots=[800];
prediction_field_snapshots=[1000 1200];
% Four absolute snapshot indices over the 3000-snapshot horizon.
error_field_snapshots=[400 800 1400 3000];
modal_atlas_n_show=4; % one conjugate side, aligned DMD/GC-DMD rows
% POD reconstruction diagnostic (reference-style panel)
n_lambda=33;
fista_opts.max_iter=5000; fista_opts.tol=1e-8; fista_opts.verbose=false;
fprintf('\n=========== GC-DMD Re=13000 cavity benchmark ===========\n');

%% 1. Load snapshots
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
lambda_grid=build_lambda_grid(train,n_lambda,'cavity13k');
[lambda_gc,~,val_gc]=optimize_operator('validation',train,val,W_train,lambda_grid,fista_opts);
[GC_raw,Pval,Sig,C_gc_final,gc_info]=estimate_gc_var1(Qmodel,gc_args{:}); %#ok<NASGU>
W_gc=build_gc_weights(C_gc_final);
A_gc=optimize_operator('solve',Qmodel,lambda_gc,W_gc,fista_opts);

%% 6. Operator diagnostics
fprintf('DMD and GC-DMD operators identified; lambda_pred=%.6g\n',lambda_gc);

% Both final operators now enter the same downstream spectral/modal analysis.
%% 7. Modal analysis: eigenvalues, frequency, growth, contribution
% Mdmd: eigenvalues, frequencies and growth rates are explicit.
[W_Mdmd,D_Mdmd]=eig(A_dmd);
lam_Mdmd=diag(D_Mdmd);
Phi_Mdmd=Ur*W_Mdmd;
b_Mdmd=W_Mdmd\Qmodel(1,:)';
freq_Mdmd=angle(lam_Mdmd)/(2*pi*dt);
growth_Mdmd=log(max(abs(lam_Mdmd),eps))/dt;
Mdmd=modal_analysis_kou(lam_Mdmd,Phi_Mdmd,b_Mdmd,n_model,dt);
Mdmd.growth_rate=growth_Mdmd(Mdmd.order);
Mdmd.frequency=freq_Mdmd(Mdmd.order); % Mgc: eigenvalues, frequencies and growth rates are explicit.
[W_Mgc,D_Mgc]=eig(A_gc);
lam_Mgc=diag(D_Mgc);
Phi_Mgc=Ur*W_Mgc;
b_Mgc=W_Mgc\Qmodel(1,:)';
freq_Mgc=angle(lam_Mgc)/(2*pi*dt);
growth_Mgc=log(max(abs(lam_Mgc),eps))/dt;
Mgc=modal_analysis_kou(lam_Mgc,Phi_Mgc,b_Mgc,n_model,dt);
Mgc.growth_rate=growth_Mgc(Mgc.order);
Mgc.frequency=freq_Mgc(Mgc.order);
fprintf('\nTop DMD modes:\n'); disp(Mdmd.table(1:min(modal_top_n,height(Mdmd.table)),:)); fprintf('Top GC-DMD modes:\n'); disp(Mgc.table(1:min(modal_top_n,height(Mgc.table)),:));

%% 8. Recursive prediction and vorticity errors for Figure 7
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
figure7_data=struct('nt',nt,'ed',ed,'eg',eg,'Xtruth',Xtruth,'Xd',Xd,'Xg',Xg);
% The helper below only prepares spatial/vorticity fields for plotting;
% time marching, reconstruction, and error evaluation remain visible here.
figure7_data=compute_cavity_figure_data(figure7_data,N,error_field_snapshots,'cavity13k');

%% 8. Paper Figure 6 (original operator/spectral drawing code)
fig6=plot_fig6(A_dmd,A_gc,Mdmd,Mgc,dt,operator_color_max,'');
%% 9. Paper Figure 7 (original recursive curve and error-field drawing code)
[fig7a,fig7b]=plot_fig7(figure7_data,n_model,pred_horizon,dt, ...
    reconstruction_field_snapshots,prediction_field_snapshots,pod_basis_modes,error_field_snapshots,time);
%% 10. Export manuscript-only results (PDF sources converted to EPS)
print(fig6,fullfile(output_dir,'fig6.eps'),'-depsc2','-painters');
print(fig7a,fullfile(output_dir,'fig7a.eps'),'-depsc2','-painters');
print(fig7b,fullfile(output_dir,'fig7b.png'),'-dpng','-r300');
fprintf('Paper Figures 6 and 7 saved in %s\n',output_dir);
