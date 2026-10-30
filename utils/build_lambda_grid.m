function lambda_grid = build_lambda_grid(data, n_lambda)
%BUILD_LAMBDA_GRID Construct a scale-aware lambda grid.
%
% The objective uses (1/(2T))*||Y-AX||_F^2.
% A natural gradient scale is max(abs(YX'/T)).

if nargin < 2
    n_lambda = 24;
end

X = data(1:end-1,:)';
Y = data(2:end,:)';
T = size(X,2);

gscale = max(abs((Y*X')/T),[],'all');
gscale = max(gscale,1e-10);

% Include zero separately; log grid spans weak to strong regularization.
lambda_grid = [0, gscale*logspace(-5,0,n_lambda-1)];

end
