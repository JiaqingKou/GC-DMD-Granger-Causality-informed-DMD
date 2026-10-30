% MAIN Reproduce manuscript Figures 1--4 for the GC-DMD VAR(1) benchmark.
% Run this file to estimate GC, identify DMD/GC-DMD/Oracle operators,
% generate the original manuscript figures, and export EPS files.

%% Step 1. Initialize environment and output directories
close all; clc; rng(42);                 % Reproducible random initialization.
here=fileparts(mfilename('fullpath'));
addpath(fullfile(here,'methods'),fullfile(here,'utils'),fullfile(here,'figures'));
output_dir=fullfile(here,'results_data');
if ~exist(output_dir,'dir'), mkdir(output_dir); end
eps_dir=fullfile(output_dir,'eps');
if ~exist(eps_dir,'dir'), mkdir(eps_dir); end

%% Step 2. Set numerical, GC, optimization, and plotting parameters
r=20; N=5000; burnin=500; sigma_process=0.02;
alpha_gc=0.05; gc_eta=1; gc_norm_mode='rowsum';
gc_intercept=true; % Include a constant in BOTH full/reduced GC regressions.
graded_density=0.30; n_lambda=32; support_tol=1e-3;
operator_clim=0.10; operator_compare_clim=0.10; gc_tail_level=0.01;
fista_opts=struct('max_iter',5000,'tol',1e-8,'verbose',false);
% taira_colors: retrieve the original manuscript's line/marker palette.
C=taira_colors();

%% Step 3. Generate the ground-truth VAR(1) operator and time series
% generate_graded_var1: build a stable sparse operator with graded couplings.
[A_true,~]=generate_graded_var1(r,0.95,graded_density,31415);
% choose_true_A_gap_threshold: obtain the nonzero-support cutoff for A_true.
true_A_tau=choose_true_A_gap_threshold(A_true);
% simulate_var1: simulate the stochastic VAR(1) trajectory with burn-in.
data=simulate_var1(A_true,N+burnin,sigma_process,1001);
data=data(burnin+1:end,:);
data=data-mean(data,1);
% Demeaning preserves the manuscript preprocessing. GC additionally fits
% an explicit intercept; demeaning alone does not enforce zero intercept.

%% Step 4. Split the sequence chronologically for identification and validation
% Training: snapshots 1:3500; validation: 3501:4250.
% The remaining snapshots are not used in this figure-generation workflow.
dataTrain=data(1:3500,:);
dataVal=data(3501:4250,:);

%% Step 5. Estimate conditional GC and identify statistically supported links
% estimate_gc_var1: fit full/reduced VAR(1) regressions (with intercept),
% calculate directional GC, and apply BH-FDR significance testing.
[GC_raw,~,Sig,C_gc,~]=estimate_gc_var1(dataTrain,'alpha',alpha_gc, ...
 'multiple_testing','bh-fdr','fdr_q',alpha_gc,'norm_mode',gc_norm_mode, ...
 'intercept',gc_intercept);
% analyze_gc_strength_plateau: analyze the normalized GC-strength tail
% for Figure 1; this is NOT the BH-FDR significance test.
G=analyze_gc_strength_plateau(GC_raw,Sig,'tail_level',gc_tail_level);
G.graph=logical(Sig)&~eye(r);
A_graph=(abs(A_true)>true_A_tau)&~eye(r);

%% Step 6. Build GC-informed weights and the regularization grid
settings=struct('alpha',alpha_gc,'norm_mode',gc_norm_mode, ...
 'eta',gc_eta,'intercept',gc_intercept);
% build_gc_weights: convert the already estimated C_gc into the penalty W;
% diagonal entries remain unpenalized (no repeated GC estimation here).
W_gc=build_gc_weights(C_gc,settings);
% build_lambda_grid: construct the data-scaled candidate lambda values.
lambda_grid=build_lambda_grid(dataTrain,n_lambda);

%% Step 7. Identify GC-DMD along the lambda path
% fit_dmd: least-squares DMD operator on the training sequence.
A_dmd_train=fit_dmd(dataTrain);
% optimize_operator('gc_path'): solve the weighted-L1 problem via FISTA
% for each lambda, evaluate validation RMSE/structural F1, select lambda,
% and refit the final GC-DMD operator using training + validation data.
Rgc=optimize_operator('gc_path',dataTrain,dataVal,W_gc,lambda_grid, ...
 A_dmd_train,A_true,support_tol,fista_opts,Sig,settings);
% choose_six_lambda_indices: select representative operators for Figures 3--4.
idx=choose_six_lambda_indices(Rgc);

%% Step 8. Identify the comparison operators (DMD and Oracle)
% fit_dmd: refit baseline DMD using training + validation snapshots.
A_dmd=fit_dmd([dataTrain;dataVal]);
A_gc=Rgc.best_A;
% Oracle weights use the known true support instead of estimated GC.
C_oracle=double(A_true~=0); C_oracle(1:r+1:end)=0;
P_oracle=1-gc_eta*C_oracle; P_oracle(1:r+1:end)=0;
% optimize_operator('validation'): choose Oracle lambda by validation RMSE.
[lambda_oracle,A_oracle]=optimize_operator('validation', ...
 dataTrain,dataVal,P_oracle,lambda_grid,fista_opts);
fprintf('GC-DMD lambda: %.5g; Oracle lambda: %.5g\n', ...
 Rgc.best_lambda,lambda_oracle);

%% Step 9. Generate the ORIGINAL manuscript plots (no styling changes)
% plot_fig1: ranked true coefficients/paired GC and two support matrices;
% also writes the Figure 1 paired-data CSV into results_data.
plot_fig1(A_true,GC_raw,G,A_graph,r,C,true_A_tau,gc_tail_level,output_dir);
% plot_fig2: four standalone operator maps (True, DMD, GC-DMD, Oracle).
plot_fig2(A_true,A_dmd,A_gc,A_oracle,r,operator_compare_clim);
% plot_fig3: six lambda-selected operators and RMSE/F1 selection curves.
plot_fig3(Rgc,idx,r,operator_clim,C);
% plot_fig4: independently ranked non-diagonal operator coefficients.
plot_fig4(Rgc,A_true,idx,support_tol);

%% Step 10. Apply manuscript formatting BEFORE exporting figures
%% This step is optional
% finalize_paper_figures: enforce the original mathematical italics,
% bold fonts, labels, and endpoint ticks across all open figures.
figs=findall(0,'Type','figure');
finalize_paper_figures(figs);

%% Step 11. Export every generated figure as a vector EPS file
% Figures are exported in creation order, with stable, descriptive names.
% print -depsc2 -painters works with older MATLAB versions; no exportgraphics
% or clim dependency is introduced. EPS export does not alter plot data.
figs=flipud(findall(0,'Type','figure'));
file_names={'Fig1a_ranked_coefficients_GC', ...
            'Fig1b_true_support', ...
            'Fig1c_GC_support', ...
            'Fig2a_true_operator', ...
            'Fig2b_DMD_operator', ...
            'Fig2c_GC_DMD_operator', ...
            'Fig2d_Oracle_operator', ...
            'Fig3_lambda_path', ...
            'Fig4_coefficient_survival'};
if numel(figs)~=numel(file_names)
    error('Expected %d figures but found %d. Check plot_fig1--4.', ...
        numel(file_names),numel(figs));
end
for k=1:numel(figs)
    fig=figs(k);
    set(fig,'PaperPositionMode','auto');
    % save as eps
    target=fullfile(eps_dir,[file_names{k} '.eps']);
    print(fig,target,'-depsc2','-painters','-r600');
    fprintf('Exported EPS (%d/%d): %s\n',k,numel(figs),target);
end
fprintf('All manuscript EPS figures saved to: %s\n',eps_dir);
