function e=one_step_relative_error(A,data)
% Relative Frobenius error of observed-state one-step predictions.
X=data(1:end-1,:)';Y=data(2:end,:)';
e=norm(Y-A*X,'fro')/max(norm(Y,'fro'),eps);
end
