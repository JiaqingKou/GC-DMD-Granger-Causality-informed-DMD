%MAIN_7 Appendix A: fixed operator, 20 supplied seeds, Figures A1--A2.
% Outcome-selected CONDITIONAL experiment, not an unbiased performance estimate.
clear; clc; close all;
root=fileparts(mfilename('fullpath'));
addpath(fullfile(root,'methods'),fullfile(root,'utils'),fullfile(root,'figures'));
output_dir=prepare_results_dir(root,7);

%% 1. Fixed experiment configuration
B=load(fullfile(root,'data','operator_bank_20.mat'));
cfg=struct('source_operator_id',13,'n_train',150,'n_val',400, ...
    'n_test',450,'burnin',500,'sigma',0.1, ...
    'n_lambda',41,'operator_color_max',0.20);
cfg.A_true=B.A_bank(:,:,cfg.source_operator_id);
cfg.seeds=[54 55 67 74 75 76 77 81 82 83 85 87 90 91 92 96 97 99 103 104];
cfg.opts=struct('max_iter',5000,'tol',1e-9,'verbose',false);

%% 2. Monte Carlo: 20 independent realizations of the same true operator
% This experiment differs from the cavity/jet tests: the known true
% operator permits direct operator-recovery metrics, not just predictions.
% All 20 realizations use the same ground-truth operator and different seeds.
A_true=cfg.A_true;
n_mc=numel(cfg.seeds);
metric_names=appendix_A_metric_names();
metrics=zeros(n_mc,3,numel(metric_names));
trials=cell(n_mc,1);
weak=abs(A_true)>0 & abs(A_true)<0.1;
zero=A_true==0;
for q=1:n_mc
    seed=cfg.seeds(q);
    data=simulate_var1(A_true,cfg.burnin+cfg.n_train+cfg.n_val+cfg.n_test,cfg.sigma,seed);
    data=data(cfg.burnin+1:end,:);
    train=data(1:cfg.n_train,:);
    val=data(cfg.n_train+1:cfg.n_train+cfg.n_val,:);
    test=data(cfg.n_train+cfg.n_val+1:end,:);

    %% 3. Standard DMD: least-squares operator
    % Fit the unregularized reference operator on training snapshots.
    % No lambda is selected for DMD; validation is used for comparisons.
    X1=train(1:end-1,:)'; X2=train(2:end,:)';
    A_dmd=X2*pinv(X1);

    %% 4. Granger causality: conditional VAR(1), intercept and BH-FDR
    % GC is learned on training data and converted to directional weights.
    % The GC-DMD diagonal is unpenalized, unlike uniform L1-DMD.
    [GC,pval,Sig,C]=estimate_gc_var1(train,'intercept',true);
    W_gc=build_gc_weights(C);  % eta=1; unpenalized diagonal

    %% 5. L1-DMD and GC-DMD: same FISTA optimizer, different penalties
    % Appendix A intentionally stores the FULL lambda paths, not just the
    % validation-selected operators: Figure A2 uses pathwise comparisons.
    % Both paths use the same candidate lambda grid and validation set;
    % their only regularization difference is the weight matrix.
    lambda=build_lambda_grid(train,cfg.n_lambda);
    W_l1=ones(size(A_true));  % Original Appendix A: diagonal penalized
    nL=numel(lambda);
    A_path=cell(nL,2);
    val_rmse=zeros(nL,2);
    test_rmse=zeros(nL,2);
    op_error=zeros(nL,2);
    active=zeros(nL,2);
    weak_error=zeros(nL,2);
    off=~eye(size(A_true,1));
    weak_off=off & weak;

    % L1-DMD: penalize every entry (including the diagonal). For each
    % lambda, retain validation/test and operator-recovery diagnostics
    % needed by the Appendix A figures; choose lambda by validation RMSE.
    for k=1:nL
        [A_candidate,info]=fit_l1_dmd(train,lambda(k),cfg.opts);
        if ~info.converged
            warning('L1-DMD FISTA did not converge: case %d, lambda %.3g.',q,lambda(k));
        end
        A_path{k,1}=A_candidate;
        val_rmse(k,1)=one_step_rmse(A_candidate,val);
        test_rmse(k,1)=one_step_rmse(A_candidate,test);
        op_error(k,1)=norm(A_candidate-A_true,'fro')/norm(A_true,'fro');
        active(k,1)=nnz(abs(A_candidate(off))>1e-3);
        weak_error(k,1)=norm(A_candidate(weak_off)-A_true(weak_off))/norm(A_true(weak_off));
    end
    [~,best_l1]=min(val_rmse(:,1));
    lambda_l1=lambda(best_l1);
    A_l1=A_path{best_l1,1};

    % GC-DMD: repeat the identical parameter sweep with GC-derived
    % weights, leaving the diagonal unpenalized. Test/true-operator errors
    % are diagnostics only and do NOT select lambda.
    for k=1:nL
        [A_candidate,info]=optimize_operator('solve',train,lambda(k),W_gc,cfg.opts);
        if ~info.converged
            warning('GC-DMD FISTA did not converge: case %d, lambda %.3g.',q,lambda(k));
        end
        A_path{k,2}=A_candidate;
        val_rmse(k,2)=one_step_rmse(A_candidate,val);
        test_rmse(k,2)=one_step_rmse(A_candidate,test);
        op_error(k,2)=norm(A_candidate-A_true,'fro')/norm(A_true,'fro');
        active(k,2)=nnz(abs(A_candidate(off))>1e-3);
        weak_error(k,2)=norm(A_candidate(weak_off)-A_true(weak_off))/norm(A_true(weak_off));
    end
    [~,best_gc]=min(val_rmse(:,2));
    lambda_gc=lambda(best_gc);
    A_gc=A_path{best_gc,2};
    best=[best_l1,best_gc];

    %% 6. Operator metrics and figure data (same selected operators)
    % Aggregate metrics use only the two validation-selected regularized
    % operators and the DMD baseline. Rcase packages the full paths solely
    % for the existing Appendix A plotting interfaces.
    As={A_dmd,A_l1,A_gc};
    mm=evaluate_operator_metrics(As,A_true,val,test,weak,zero);
    metrics(q,:,:)=reshape(mm,[1 3 numel(metric_names)]);
    % The original Appendix A plotting functions expect these fields.
    Rcase=struct('DMD',A_dmd,'A',{A_path},'best',best, ...
        'lambda',lambda,'val',val_rmse,'test',test_rmse, ...
        'operator_error',op_error,'active',active,'weak_error',weak_error, ...
        'GC',GC,'p',pval,'Sig',Sig,'C',C,'W1',W_l1,'Wgc',W_gc);
    trials{q}=struct('seed',seed,'A_true',A_true,'R',Rcase);
    fprintf('Case %2d/%d, seed=%d: L1/GC operator error %.6g/%.6g; lambda %.3g/%.3g\n', ...
        q,n_mc,seed,mm(2,1),mm(3,1),lambda_l1,lambda_gc);
end
mean_metrics=squeeze(mean(metrics,1));
summary_table=array2table(mean_metrics,'VariableNames',metric_names, ...
    'RowNames',{'DMD','L1-DMD','GC-DMD'});
disp('Appendix A: mean metrics across 20 realizations:');
disp(summary_table);
R=struct('metrics',metrics,'trials',{trials});

%% 7. Representative case (same selection rule as the original project)
% The displayed example maximizes relative GC-DMD operator-error
% improvement over L1-DMD; the 20-run averages above use ALL realizations.
errors=R.metrics(:,:,1);
improvement=(errors(:,2)-errors(:,3))./max(errors(:,2),eps);
[~,idx]=max(improvement);
t=R.trials{idx};
fprintf('Appendix A: representative seed %d (selected by relative operator-error reduction).\n',t.seed);

%% 8. Original manuscript figures and EPS names
figA1=plot_figA1(t.R,t.A_true,cfg.operator_color_max);
for j=1:4
    print(figA1(j),fullfile(output_dir,sprintf('figA1%c.eps',char('a'+j-1))),'-depsc2','-painters');
end
figA2=plot_figA2(t.R);
print(figA2,fullfile(output_dir,'figA2.eps'),'-depsc2','-painters');
fprintf('Appendix A figures saved to %s\n',output_dir);
