function R=compute_cavity_pod(X,Xc,mu,N,time,Re,cfg)
%COMPUTE_CAVITY_POD Figure 5 postprocessing; input snapshots already centered.
[U,s,V]=compute_svd(Xc);
d=s.^2; energy=d/sum(d); cumulative=cumsum(energy);
r=find(cumulative>=cfg.pod_energy_target,1);
q1=s(1)*V(:,1);
[~,anchor]=max(abs(U(:,1)));
if U(anchor,1)<0, q1=-q1; end
scale=max(abs(q1)); Qnorm=q1/scale;
idx=cfg.snapshot_indices;
assert(all(idx>=1 & idx<=cfg.n_snapshots & idx==fix(idx)),'Invalid snapshot indices.');
fields=X(:,idx);
if strcmp(cfg.field_source,'pod')
    U=bsxfun(@rdivide,Xc*V(:,1:r),sqrt(d(1:r))');
    fields=bsxfun(@plus,mu,U*(U'*Xc(:,idx)));
elseif ~strcmp(cfg.field_source,'truth')
    error('field_source must be truth or pod.');
end
% Same -del2 as the original toolbox, evaluated as tensor products to avoid
% building a dense ((N+1)^2 x (N+1)^2) Laplacian solely for three fields.
Grid=CavityGridOperators(N); D2=Grid.D^2;
[~,ord]=sort(Grid.x); x=Grid.x(ord);
xi=linspace(-1,1,cfg.interpolation_points); [xxi,yyi]=meshgrid(xi,xi);
W=cell(1,numel(idx));
for j=1:numel(idx)
    Psi=reshape(fields(:,j),N+1,N+1);
    vort=-(D2*Psi+Psi*D2');
    % Ascending coordinates prevent release-dependent interp2 handling of descending grids.
    W{j}=interp2(x,x,vort(ord,ord),xxi,yyi,'spline');
    assert(all(isfinite(W{j}(:))),'Interpolation failed for Re=%d snapshot=%d.',Re,idx(j));
end
R=struct('Re',Re,'N',N,'rank',r,'retained_energy',cumulative(r), ...
    'first_energy',energy(1),'energy',energy,'cumulative_energy',cumulative, ...
    'q1',q1,'q1_normalized',Qnorm,'normalization_scale',scale, ...
    'time',time,'snapshot_indices',idx,'field_source',cfg.field_source, ...
    'xi',xi,'vorticity',{W});
fprintf('Re=%d: r=%d, retained=%.8f%%, q1 maxabs=%.6g\n',Re,r,100*cumulative(r),scale);
end
