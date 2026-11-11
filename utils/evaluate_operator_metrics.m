function mm=evaluate_operator_metrics(As,A_true,val,test,weak,zero)
%EVALUATE_OPERATOR_METRICS Common metrics for a known reference operator.
if nargin<5, weak=abs(A_true)>0 & abs(A_true)<0.1; end
if nargin<6, zero=A_true==0; end
mm=zeros(numel(As),9);
for j=1:numel(As)
    A=As{j}; E=A-A_true;
    mm(j,:)=[norm(E,'fro')/norm(A_true,'fro'), ...
        one_step_rmse(A,val),one_step_rmse(A,test), ...
        norm(E(weak))/norm(A_true(weak)),sqrt(mean(E(zero).^2)), ...
        sqrt(mean(diag(E).^2)),max(abs(eig(A))), ...
        one_step_relative_error(A,val),one_step_relative_error(A,test)];
end
end
