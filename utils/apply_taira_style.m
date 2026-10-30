function apply_taira_style(obj)
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
