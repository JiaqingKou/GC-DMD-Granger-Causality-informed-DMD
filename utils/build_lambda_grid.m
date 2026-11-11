function lambda_grid = build_lambda_grid(data, n_lambda, case_name)
%BUILD_LAMBDA_GRID Construct a scale-aware lambda grid.
%
% The objective uses (1/(2T))*||Y-AX||_F^2.
% A natural gradient scale is max(abs(YX'/T)).

if nargin >= 3 && (ischar(case_name) || isstring(case_name)) && strcmpi(char(case_name),'jet')
    lambda_grid = jet_lambda_grid_impl(data,n_lambda,[-10 1]);
    return;
end
if nargin >= 3 && ismember(lower(string(case_name)),["cavity13k","cavity19k","cavity30k"])
    lambda_grid = cavity13k_lambda_impl(data);
    return;
end
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

function lambda_grid = cavity13k_lambda_impl(data, n_lambda)
%BUILD_LAMBDA_GRID Construct a nested, scale-aware lambda grid.
%
% The objective uses (1/(2T))*||Y-AX||_F^2.
% A natural gradient scale is max(abs(YX'/T)).
%
% The logarithmic grid uses a fixed number of points per decade. Therefore,
% extending lambda_exp_min only ADDS weaker candidates and does not move any
% existing candidates.  The input n_lambda is retained only for backward
% compatibility with the existing project interface.

if nargin < 2
    n_lambda = []; %#ok<NASGU>
else
    n_lambda = n_lambda; %#ok<NASGU>
end

X = data(1:end-1,:)';
Y = data(2:end,:)';
T = size(X,2);

gscale = max(abs((Y*X')/T),[],'all');
gscale = max(gscale,1e-10);

% Nested lambda-grid settings.
lambda_exp_min = -12;
lambda_exp_max = 0;
points_per_decade = 8;

exponents = lambda_exp_min : 1/points_per_decade : lambda_exp_max;
lambda_grid = [0, gscale * 10.^exponents];

fprintf('Lambda grid: zero + [%.3e, %.3e], %d nonzero candidates (%d/decade)\n', ...
    lambda_grid(2), lambda_grid(end), numel(lambda_grid)-1, points_per_decade);

end

function lambda_grid = jet_lambda_grid_impl(data, n_lambda, exp_range)
%BUILD_LAMBDA_GRID Scale-aware wide lambda grid for GC-DMD tuning.
% Objective: (1/(2T))*||Y-AX||_F^2 + lambda*penalty.
% Default spans 11 decades around the data gradient scale so that both the
% under-regularized plateau and over-regularized branch are visible.
if nargin < 2 || isempty(n_lambda), n_lambda = 67; end
if nargin < 3 || isempty(exp_range), exp_range = [-10 1]; end
X=data(1:end-1,:)'; Y=data(2:end,:)'; T=size(X,2);
gscale=max(abs((Y*X')/T),[],'all'); gscale=max(gscale,1e-14);
lambda_grid=[0, gscale*logspace(exp_range(1),exp_range(2),n_lambda-1)];
fprintf('Lambda grid: zero + [%.3e, %.3e], %d nonzero candidates\n',lambda_grid(2),lambda_grid(end),numel(lambda_grid)-1);
end
