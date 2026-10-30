function A = fit_dmd(data)
%FIT_DMD Ordinary least-squares DMD / VAR(1) operator.
%
% data : N x r
%
% x_{k+1} = A x_k

X1 = data(1:end-1,:)';
X2 = data(2:end,:)';

A = X2 * pinv(X1);

end
