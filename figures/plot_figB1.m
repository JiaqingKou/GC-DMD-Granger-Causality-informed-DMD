function fig=plot_figB1(G,lambda_grid,val_gc,val_dmd,lambda_pred,Ad,Ag,cmax,gmax,outdir,validation_ylim)
% One canvas: GC | lambda selection; DMD operator | GC-DMD operator.
if nargin<11,validation_ylim=[0 .012];end
C=taira_colors(); fig=figure('Name','GC selection and operators','Color','w','Position',[80 60 1100 1100]);
pp={[.10 .60 .30 .30],[.59 .60 .30 .30],[.10 .12 .30 .30],[.59 .12 .30 .30]};
axs=gobjects(1,4);
for k=1:4
    ax=axes(fig,'Position',pp{k}); axs(k)=ax;
    if k==2
        pos=lambda_grid>0;
        semilogx(ax,lambda_grid(pos),val_gc(pos),'-o','Color',C.red,'LineWidth',1.5,'MarkerSize',3); hold(ax,'on');
        yline(ax,val_dmd,'--','Color',C.gray,'LineWidth',1.1);
        ylim(ax,validation_ylim); yticks(ax,linspace(validation_ylim(1),validation_ylim(2),7));
        if lambda_pred>0
            xline(ax,lambda_pred,'--','Color',C.blue,'HandleVisibility','off');
            yl=ylim(ax);
            text(ax,lambda_pred,0.006,'\lambda_{pred}', ...
                'Interpreter','tex','Rotation',90,'VerticalAlignment','middle', ...
                'HorizontalAlignment','center', ...
                'Color',C.blue,'FontSize',18,'Clipping','on');
        end
        lp=lambda_grid(pos);
        lo=floor(log10(min(lp))); hi=ceil(log10(max(lp)));
        xlim(ax,[10^lo 10^hi]);
        exponents=unique([lo, lo+3, lo+6, lo+9, hi]);
        exponents=exponents(exponents>=lo & exponents<=hi);
        xticks(ax,unique([10^lo 10.^exponents 10^hi])); % include both axis endpoints
        ax.XTickLabelRotation=0;
        xlabel(ax,'\lambda'); ylabel(ax,'RMSE_{val}(\itA\rm)');
        legend(ax,'GC-DMD','DMD (\lambda=0)','Location','best');
    else
        if k==1
            imagesc(ax,G); caxis(ax,[0 gc_display_limit(G,gmax)]); colormap(ax,gc_white_red(256));
            xlabel(ax,'Cause \itj'); ylabel(ax,'Effect \iti'); tag='GC';
        else
            if k==3, A=Ad; else,A=Ag;end
            imagesc(ax,A); caxis(ax,[-cmax cmax]); colormap(ax,taira_diverging(256));
            xlabel(ax,'Source \itj'); ylabel(ax,'Target \iti'); tag='A_{ij}';
        end
        axis(ax,'image'); ax.Tag='PaperMatrixAxes'; paper_matrix_ticks(ax,size(G,1));
        cb=colorbar(ax); cb.Position=[pp{k}(1)+.315 pp{k}(2) .012 .30]; cb.Ticks=linspace(ax.CLim(1),ax.CLim(2),5);

    end
end
apply_taira_style(fig);
for k=1:4
    axs(k).PositionConstraint='innerposition'; axs(k).Position=pp{k}; pbaspect(axs(k),[1 1 1]);
    % A second xlabel line reserves space for the panel letter.
    axis_label=axs(k).XLabel.String;
    axs(k).XLabel.String={axis_label,sprintf('(%c)','a'+k-1)};
end
legs=findall(fig,'Type','legend'); for k=1:numel(legs),legs(k).FontSize=17;end
cbs=findall(fig,'Type','colorbar'); for k=1:numel(cbs),cbs(k).FontSize=17;end
% Defer graphics rendering; no forced refresh.
for k=[1 3 4]
    axs(k).Position=pp{k};
end
set(fig,'PaperUnits','inches','PaperPositionMode','manual','PaperPosition',[0 0 8 7.2],'PaperSize',[8 7.2]);
% Display only; no output folder required.
% Automatic EPS export disabled; display first.
% No automatic PDF export.
end


function a=gc_display_limit(G,requested)
% Display only: omit self terms; use 99th percentile of positive entries.
if nargin>1 && ~isempty(requested), a=requested; else
    v=G(~eye(size(G))); v=v(isfinite(v)&v>0);
    if isempty(v), a=1; else, a=prctile(v,99); end
end
if ~isscalar(a)||~isfinite(a)||a<=0, error('GC color maximum must be positive.'); end
end


function cmap = gc_white_red(n)
%GC_WHITE_RED White at zero; increasingly saturated red for positive GC.
if nargin < 1 || isempty(n), n = 256; end
n = max(2,round(n));
t = linspace(0,1,n)';
% Match the operator-style white background and red active couplings.
anchor = [1 1 1; 0.98 0.89 0.88; 0.91 0.53 0.51; 0.70 0.11 0.12];
cmap = interp1(linspace(0,1,size(anchor,1)),anchor,t,'linear');
end


function paper_matrix_ticks(ax,n)
% At most six integer ticks, including both endpoints.
t=unique(round(linspace(1,n,min(n,6))));
set(ax,'XTick',t,'YTick',t,'XTickLabelRotation',0,'YTickLabelRotation',0);
end


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
    % Explicitly show both axis limits as major tick values.
    % Do not modify image axes whose integer ticks are handled separately.
    try
        if ~strcmp(ax.Tag,'PaperMatrixAxes') && strcmp(ax.XScale,'linear') && ~isempty(ax.XTick)
            xl=ax.XLim; ax.XTick=unique([xl(1),ax.XTick(ax.XTick>xl(1) & ax.XTick<xl(2)),xl(2)]);
        end
        if ~strcmp(ax.Tag,'PaperMatrixAxes') && strcmp(ax.YScale,'linear') && ~isempty(ax.YTick)
            yl=ax.YLim; ax.YTick=unique([yl(1),ax.YTick(ax.YTick>yl(1) & ax.YTick<yl(2)),yl(2)]);
        end
    catch
    end
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
if isappdata(fig,'PaperPanelLayout') && getappdata(fig,'PaperPanelLayout')
    align_paper_panel(fig);
end
% Styling does not force all open figures to refresh.
end


function align_paper_panel(fig)
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


function C = taira_colors()
%TAIRA_COLORS Consistent, colorblind-conscious palette for publication figures.
C.red = [0.8500 0.3250 0.0980];
C.blue   = [0.0000 0.4470 0.7410];
C.orange = [0.8500 0.3250 0.0980];
C.yellow = [0.9290 0.6940 0.1250];
C.purple = [0.4940 0.1840 0.5560];
C.green  = [0.4660 0.6740 0.1880];
C.cyan   = [0.3010 0.7450 0.9330];
C.gray   = [0.35 0.35 0.35];
C.lightgray = [0.72 0.72 0.72];
end
