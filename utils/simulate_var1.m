function data = simulate_var1(A, N, sigma_process, seed)
%SIMULATE_VAR1 Simulate x_{k+1}=A x_k + epsilon_k.
%
% Output:
%   data : N x r

if nargin < 3
    sigma_process = 0.02;
end
if nargin < 4
    seed = 1;
end

rng(seed);

r = size(A,1);
X = zeros(r,N);
X(:,1) = randn(r,1);

for k = 1:N-1
    eps_k = sigma_process * randn(r,1);
    X(:,k+1) = A*X(:,k) + eps_k;
end

data = X';

end
