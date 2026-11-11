function apply_taira_style(obj,case_name)
if nargin>=2 && strcmpi(case_name,'cavity19k')
    cavity19k_style_impl(obj);
    return;
end
if nargin>=2 && strcmpi(case_name,'cavity13k')
    cavity13k_style_impl(obj);
    return;
end
%APPLY_TAIRA_STYLE Lightweight Taira-like publication styling.
% Times New Roman, larger/bold labels, restrained grid, no graphics listeners.
if nargin<1 || isempty(obj), obj=gcf; end
if isa(obj,'matlab.ui.Figure'), fig=obj; else, fig=ancestor(obj,'figure'); end
if isempty(fig), fig=gcf; end
set(fig,'Color','w','Renderer','opengl');
axs=findall(fig,'Type','axes');
for k=1:numel(axs)
    ax=axs(k);
    set(ax,'FontName','Times New Roman','FontSize',20,'FontWeight','bold', ...
        'LineWidth',1.0,'TickDir','out','Box','on','Layer','top', ...
        'XMinorTick','off','YMinorTick','off');
    grid(ax,'off'); set(ax,'XGrid','off','YGrid','off','XMinorGrid','off','YMinorGrid','off');
    ax.TickLabelInterpreter='tex';
    ax.XLabel.Interpreter='tex'; ax.YLabel.Interpreter='tex';
    ax.XLabel.FontName='Times New Roman'; ax.YLabel.FontName='Times New Roman';
    ax.XLabel.FontSize=22; ax.YLabel.FontSize=22;
    ax.XLabel.FontWeight='bold'; ax.YLabel.FontWeight='bold';
    ax.Title.Interpreter='tex'; ax.Title.FontName='Times New Roman';
    ax.Title.String=''; ax.Title.FontSize=16; ax.Title.FontWeight='bold';
end
legs=findall(fig,'Type','legend');
for k=1:numel(legs)
    set(legs(k),'FontName','Times New Roman','FontSize',16,'FontWeight','bold', ...
        'Box','off','Interpreter','tex');
end
cbs=findall(fig,'Type','colorbar');
for k=1:numel(cbs)
    set(cbs(k),'FontName','Times New Roman','FontSize',16,'FontWeight','bold','TickDirection','out');
    try
        cbs(k).Label.Interpreter='tex'; cbs(k).Label.FontName='Times New Roman';
        cbs(k).Label.FontSize=20; cbs(k).Label.FontWeight='bold';
    catch
    end
end
drawnow limitrate;
end

function cavity13k_style_impl(obj)
%APPLY_TAIRA_STYLE Lightweight Taira-like publication styling.
% Times New Roman, larger/bold labels, restrained grid, no graphics listeners.
if nargin<1 || isempty(obj), obj=gcf; end
if isa(obj,'matlab.ui.Figure'), fig=obj; else, fig=ancestor(obj,'figure'); end
if isempty(fig), fig=gcf; end
set(fig,'Color','w','Renderer','opengl');
axs=findall(fig,'Type','axes');
for k=1:numel(axs)
    ax=axs(k);
    if strcmp(ax.Tag,'PaperRowLabels'), continue; end
    set(ax,'FontName','Times New Roman','FontSize',13,'FontWeight','bold', ...
        'LineWidth',1.0,'TickDir','out','Box','on','Layer','top', ...
        'XMinorTick','off','YMinorTick','off');
    grid(ax,'off');
    set(ax,'XGrid','off','YGrid','off','ZGrid','off', ...
        'XMinorGrid','off','YMinorGrid','off','ZMinorGrid','off');
    ax.TickLabelInterpreter='tex';
    ax.XLabel.Interpreter='tex'; ax.YLabel.Interpreter='tex';
    ax.XLabel.FontName='Times New Roman'; ax.YLabel.FontName='Times New Roman';
    ax.XLabel.FontSize=15; ax.YLabel.FontSize=15;
    ax.XLabel.FontWeight='bold'; ax.YLabel.FontWeight='bold';
    ax.Title.String='';
    ax.Title.Interpreter='tex'; ax.Title.FontName='Times New Roman';
    ax.Title.FontSize=16; ax.Title.FontWeight='bold';
end
legs=findall(fig,'Type','legend');
for k=1:numel(legs)
    set(legs(k),'FontName','Times New Roman','FontSize',12,'FontWeight','bold', ...
        'Box','off','Interpreter','tex');
end
cbs=findall(fig,'Type','colorbar');
for k=1:numel(cbs)
    set(cbs(k),'FontName','Times New Roman','FontSize',12,'FontWeight','bold','TickDirection','out');
    try
        cbs(k).Label.Interpreter='tex'; cbs(k).Label.FontName='Times New Roman';
        cbs(k).Label.FontSize=14; cbs(k).Label.FontWeight='bold';
    catch
    end
end
layouts=findall(fig,'Type','tiledlayout');
for k=1:numel(layouts), layouts(k).Title.String=''; end
sg=findall(fig,'Tag','suptitle'); delete(sg);
% Enlarged fonts only for the requested operator/spectral and story figures.
paper_names={'DMD operator','GC-DMD operator','Ritz eigenvalues', ...
    'Frequency contribution spectrum','Combined operator and spectral paper figure', ...
    'Recursive reconstruction and prediction','Signed vorticity error slices'};
if any(strcmp(fig.Name,paper_names))
    for ia=1:numel(axs)
        if strcmp(axs(ia).Tag,'PaperRowLabels'),continue;end
        set(axs(ia),'FontSize',20,'XTickLabelRotation',0);
        set([axs(ia).XLabel axs(ia).YLabel],'FontSize',22);
    end
    set(legs,'FontSize',17);
    for ic=1:numel(cbs),cbs(ic).FontSize=17;cbs(ic).Label.FontSize=19;end
end
% Balanced paper typography: larger than the original 13/15 pt profile.
if strcmp(fig.Name,'Signed vorticity error slices')
    for ia=1:numel(axs)
        if strcmp(axs(ia).Tag,'PaperRowLabels'),continue;end
        set(axs(ia),'FontSize',20,'XTick',[-1 0 1],'YTick',[-1 0 1]);
        set([axs(ia).XLabel axs(ia).YLabel],'FontSize',22);
        tt=findall(axs(ia),'Type','text');
        for it=1:numel(tt)
            if (startsWith(string(tt(it).String),'k =') || startsWith(string(tt(it).String),'\it k')),tt(it).FontSize=18;end
        end
    end
    for ic=1:numel(cbs),cbs(ic).FontSize=18;cbs(ic).Label.FontSize=24;end
end
if isappdata(fig,'PaperPanelLayout') && getappdata(fig,'PaperPanelLayout')
    align_paper_panel(fig);
end
% Defer rendering to MATLAB; styling must not force all open figures to redraw.
end

function cavity19k_style_impl(obj)
%APPLY_TAIRA_STYLE Lightweight Taira-like publication styling.
% Times New Roman, larger/bold labels, restrained grid, no graphics listeners.
if nargin<1 || isempty(obj), obj=gcf; end
if isa(obj,'matlab.ui.Figure'), fig=obj; else, fig=ancestor(obj,'figure'); end
if isempty(fig), fig=gcf; end
set(fig,'Color','w','Renderer','opengl');
axs=findall(fig,'Type','axes');
for k=1:numel(axs)
    ax=axs(k);
    if strcmp(ax.Tag,'PaperRowLabels'), continue; end
    set(ax,'FontName','Times New Roman','FontSize',20,'FontWeight','bold', ...
        'LineWidth',1.0,'TickDir','out','Box','on','Layer','top', ...
        'XMinorTick','off','YMinorTick','off');
    grid(ax,'off');
    set(ax,'XGrid','off','YGrid','off','ZGrid','off', ...
        'XMinorGrid','off','YMinorGrid','off','ZMinorGrid','off');
    ax.TickLabelInterpreter='tex';
    ax.XLabel.Interpreter='tex'; ax.YLabel.Interpreter='tex';
    ax.XLabel.FontName='Times New Roman'; ax.YLabel.FontName='Times New Roman';
    ax.XLabel.FontSize=22; ax.YLabel.FontSize=22;
    ax.XLabel.FontWeight='bold'; ax.YLabel.FontWeight='bold';
    ax.XTickLabelRotation=0; ax.YTickLabelRotation=0;
    tt=findall(ax,'Type','text');
    for it=1:numel(tt)
        if ~isequal(tt(it),ax.XLabel) && ~isequal(tt(it),ax.YLabel) && ~isequal(tt(it),ax.Title)
            set(tt(it),'FontName','Times New Roman','FontSize',18,'FontWeight','bold');
        end
    end
    ax.Title.String='';
    ax.Title.Interpreter='tex'; ax.Title.FontName='Times New Roman';
    ax.Title.FontSize=16; ax.Title.FontWeight='bold';
end
legs=findall(fig,'Type','legend');
for k=1:numel(legs)
    set(legs(k),'FontName','Times New Roman','FontSize',17,'FontWeight','bold', ...
        'Box','off','Interpreter','tex');
end
cbs=findall(fig,'Type','colorbar');
for k=1:numel(cbs)
    set(cbs(k),'FontName','Times New Roman','FontSize',17,'FontWeight','bold','TickDirection','out');
    try
        cbs(k).Label.Interpreter='tex'; cbs(k).Label.FontName='Times New Roman';
        cbs(k).Label.FontSize=19; cbs(k).Label.FontWeight='bold';
    catch
    end
end
layouts=findall(fig,'Type','tiledlayout');
for k=1:numel(layouts), layouts(k).Title.String=''; end
sg=findall(fig,'Tag','suptitle'); delete(sg);
% Slightly larger legends for rank-sensitivity paper panels.
if startsWith(fig.Name,'Rank robustness')
    set(legs,'FontSize',19);
end
if isappdata(fig,'PaperPanelLayout') && getappdata(fig,'PaperPanelLayout')
    align_paper_panel(fig,'cavity19k');
end
% Defer rendering to MATLAB to avoid blocking all open figures.
end
