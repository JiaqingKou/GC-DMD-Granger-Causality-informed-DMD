function [A, info] = generate_graded_var1(r,target_rho,density,seed)
%GENERATE_GRADED_VAR1 Stable VAR(1) with moderately dense, continuously graded couplings.
%
% This case is intended for the |A_ij|-versus-GC threshold diagnostic.
% Unlike generate_sparse_var1, off-diagonal coefficients do not have a
% sharp nonzero/zero magnitude gap.  Their magnitudes are continuously
% distributed from weak to strong values.

if nargin < 1, r = 20; end
if nargin < 2, target_rho = 0.95; end
if nargin < 3, density = 0.30; end
if nargin < 4, seed = 31415; end
if density <= 0 || density > 1, error('density must be in (0,1].'); end
rng(seed);

A = zeros(r);
% Heterogeneous self-memory, same order as cross-couplings.
diagv = 0.08 + 0.55*rand(r,1);
A(1:r+1:end) = diagv;

off = find(~eye(r));
nOff = numel(off);
nEdge = max(1,round(density*nOff));
sel = off(randperm(nOff,nEdge));

% Log-uniform magnitudes give a smooth weak-to-strong hierarchy.
mag_min = 0.008;
mag_max = 0.35;
u = rand(nEdge,1);
mag = mag_min .* (mag_max/mag_min).^u;
signv = 2*(rand(nEdge,1)>0.5)-1;
A(sel) = signv .* mag;

% Stabilize globally while preserving the relative coefficient hierarchy.
rho0 = max(abs(eig(A)));
if rho0 > 0
    A = A * (target_rho/rho0);
end

info.density_requested = density;
info.offdiag_edges = nEdge;
info.offdiag_density = nEdge/nOff;
info.rho_before_scaling = rho0;
info.rho_after_scaling = max(abs(eig(A)));
info.description = 'graded dense VAR with log-uniform off-diagonal magnitudes';
end
