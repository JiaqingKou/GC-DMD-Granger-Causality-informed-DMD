function align_paper_panel(fig,case_name)
if nargin>=2 && strcmpi(case_name,'cavity19k')
    align_cavity19k_panel(fig);
    return;
end
% Identical canvas and square plot box; reserve colorbar space in all panels.
setappdata(fig,'PaperPanelLayout',true);
fig.Units='pixels'; pos=fig.Position; pos(3:4)=[850 850]; fig.Position=pos;
axs=findall(fig,'Type','axes');
for k=1:numel(axs)
    if strcmp(axs(k).Tag,'PaperRowLabels'), continue; end
    axs(k).Units='normalized';
    axs(k).PositionConstraint='innerposition';
    axs(k).Position=[0.20 0.20 0.60 0.60];
    pbaspect(axs(k),[1 1 1]);
end
cbs=findall(fig,'Type','colorbar');
for k=1:numel(cbs)
    cbs(k).Units='normalized'; cbs(k).Position=[0.84 0.20 0.025 0.60];
end
% Reapply after colorbar positioning to prevent automatic axes shrinkage.
for k=1:numel(axs)
    if ~strcmp(axs(k).Tag,'PaperRowLabels'), axs(k).Position=[0.20 0.20 0.60 0.60]; end
end
fig.PaperUnits='inches'; fig.PaperPositionMode='manual';
fig.PaperPosition=[0 0 8.5 8.5]; fig.PaperSize=[8.5 8.5];
end

function align_cavity19k_panel(fig)
% Identical canvas and square plot box; reserve colorbar space in all panels.
setappdata(fig,'PaperPanelLayout',true);
fig.Units='pixels'; pos=fig.Position; pos(3:4)=[600 600]; fig.Position=pos;
axs=findall(fig,'Type','axes');
for k=1:numel(axs)
    if strcmp(axs(k).Tag,'PaperRowLabels'), continue; end
    axs(k).Units='normalized';
    axs(k).PositionConstraint='innerposition';
    axs(k).Position=[0.16 0.16 0.68 0.68];
    pbaspect(axs(k),[1 1 1]);
end
cbs=findall(fig,'Type','colorbar');
for k=1:numel(cbs)
    cbs(k).Units='normalized'; cbs(k).Position=[0.865 0.16 0.028 0.68];
end
% Reapply after colorbar positioning to prevent automatic axes shrinkage.
for k=1:numel(axs)
    if ~strcmp(axs(k).Tag,'PaperRowLabels'), axs(k).Position=[0.16 0.16 0.68 0.68]; end
end
fig.PaperUnits='inches'; fig.PaperPositionMode='manual';
fig.PaperPosition=[0 0 6 6]; fig.PaperSize=[6 6];
end
