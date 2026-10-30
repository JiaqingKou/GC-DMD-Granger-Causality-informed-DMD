function rmse = one_step_rmse(A, data)
%ONE_STEP_RMSE One-step RMSE over all states and all available transitions.

if size(data,1) < 2
    rmse = NaN;
    return;
end

X1 = data(1:end-1,:)';
X2 = data(2:end,:)';

E = X2 - A*X1;
rmse = sqrt(mean(E(:).^2));

end
