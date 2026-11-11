function draw_reference_flow(ax,x,y,W,limits)
% Signed reference display: retain the pale-blue cavity background.
% Scalar legacy limits are replaced by a signed percentile range.
if nargin<5 || isempty(limits) || isscalar(limits)
    v=real(W(isfinite(W))); limits=prctile(v,[1 99]);
end
if numel(limits)~=2 || any(~isfinite(limits)) || limits(2)<=limits(1)
    limits=[min(real(W(:))) max(real(W(:)))];
    if limits(2)<=limits(1), limits=limits(1)+[-1 1]; end
end
Wp=max(min(real(W),limits(2)),limits(1));
contourf(ax,x,y,Wp,linspace(limits(1),limits(2),101),'LineStyle','none');
axis(ax,'square'); xlim(ax,[-1 1]); ylim(ax,[-1 1]);
caxis(ax,limits); colormap(ax,VorticityColormap(64));
set(ax,'XTick',[-1 0 1],'YTick',[-1 0 1]); box(ax,'on'); grid(ax,'off');
xlabel(ax,'\itx'); ylabel(ax,'\ity');
end
