function plot_fig5(Results,cfg)
%PLOT_FIG5 Manuscript Figure 5: original independent plotting windows.
% The two original plotting implementations are preserved below as local functions.
plot_pod1_comparison(Results,cfg);
plot_snapshot_comparison(Results,cfg);
end

function fig=plot_pod1_comparison(Results,cfg)
% Overlay the three cases, as used in the manuscript.
fig=figure('Color','w','Name','Normalized first POD coefficients','NumberTitle','off', ...
    'Position',cfg.coefficient_position);
ax=axes(fig,'Units','normalized','Position',[0.12 0.25 0.84 0.61]);hold(ax,'on');
cols=[0.12 0.32 0.70;0.85 0.32 0.16;0 0 0];styles={'-','--','-.'};
h=gobjects(1,numel(Results));labels=cell(1,numel(Results));
for j=1:numel(Results)
    R=Results{j};
    h(j)=plot(ax,1:cfg.n_snapshots,R.q1_normalized,'Color',cols(j,:), ...
        'LineStyle',styles{j},'LineWidth',cfg.line_width);
    labels{j}=['\bf\it Re\rm\bf = ' num2str(R.Re)];
end
xlim(ax,[1 cfg.n_snapshots]);ylim(ax,[-1.2 1.2]);
for j=1:numel(cfg.snapshot_indices)
    k=cfg.snapshot_indices(j);
    xline(ax,k,'--','Color',[0.95 0.70 0.10],'LineWidth',1.6,'HandleVisibility','off');
    ha='center';if j==1,ha='left';elseif j==numel(cfg.snapshot_indices),ha='right';end
    text(ax,k,1.08,sprintf('\\fontname{Times New Roman}\\bf\\it k\\rm\\bf = %d',k),'HorizontalAlignment',ha,'VerticalAlignment','bottom', ...
        'FontName',cfg.font_name,'FontSize',cfg.annotation_size,'FontWeight','bold', ...
        'BackgroundColor','w','Margin',2,'HandleVisibility','off','Interpreter','tex');
end
set(ax,'XTick',[1 100 200 300 400 500 600 700 800],'YTick',-1.2:0.6:1.2);
ylabel(ax,'\fontname{Times New Roman}\bf\it q\rm\bf_1 / max|\it q\rm\bf_1|','Interpreter','tex');
xlabel(ax,['\bf\it k\rm\bf (\Delta\it t\rm\bf = ' num2str(cfg.nondimensional_dt,'%.1f') ')'], ...
    'Interpreter','tex','FontName',cfg.font_name,'FontWeight','bold');
style_cavity_axes(ax,cfg);set(ax,'XTickLabelRotation',0);
lg=legend(ax,h,labels,'Location','northoutside','Orientation','horizontal','AutoUpdate','off');
set(lg,'FontName',cfg.font_name,'FontSize',cfg.legend_size,'FontWeight','bold','FontAngle','normal','Interpreter','tex','Box','off');
% Restore plot box after legend layout; reserve bottom room for the two-line label.
ax.Position=[0.12 0.25 0.84 0.61];
export_cavity_figure(fig,'fig5a',cfg);
end

function fig=plot_snapshot_comparison(Results,cfg)
% One colorbar per Reynolds number. Three times in one row share physical limits.
nr=numel(Results); nc=numel(cfg.snapshot_indices);
assert(nr==3 && nc==3,'This figure layout expects three cases and three snapshots.');
limits=zeros(nr,1);
for k=1:nr
    if isempty(cfg.field_color_limit)
        values=[];
        for j=1:nc, values=[values;abs(Results{k}.vorticity{j}(:))]; end %#ok<AGROW>
        values=sort(values(isfinite(values)));
        p=cfg.field_color_percentile;
        assert(p>=0 && p<=100,'Percentile must lie in [0,100].');
        pos=1+(numel(values)-1)*p/100;
        lo=floor(pos); hi=ceil(pos); a=values(lo)+(pos-lo)*(values(hi)-values(lo));
        if a<=0, a=max(values); end
        if a<=0, a=1; end
    elseif isscalar(cfg.field_color_limit)
        a=cfg.field_color_limit;
    else
        assert(numel(cfg.field_color_limit)==nr,'Specify one limit or three Reynolds-specific limits.');
        a=cfg.field_color_limit(k);
    end
    assert(isscalar(a) && isfinite(a) && a>0,'Invalid vorticity color limit.');
    limits(k)=a;
end
if strcmp(cfg.field_colormap,'original')
    cmap=VorticityColormap(64);
else
    cmap=cavity_smooth_redblue(257);
end
% Fixed paper canvas/fonts; leave the top toolbar accessible.
fp=cfg.field_position; screen=get(groot,'ScreenSize');
fw=fp(3); left=110; right=60; top=90; bottom=115;
gapx=105; gapy=215;
% Derive height from the actual three rows, not an oversized empty canvas.
side=(fw-left-right-(nc-1)*gapx)/nc;
assert(side>0,'Field canvas is too small.');
fh=top+nr*side+(nr-1)*gapy+bottom;
fp(4)=fh; fp(2)=screen(2)+screen(4)-fh-110;
fig=figure('Color','w','Name','Cavity snapshots','NumberTitle','off', ...
    'Units','pixels','WindowStyle','normal','MenuBar','figure','ToolBar','figure','Position',fp);
set(fig,'Position',fp);
setappdata(fig,'CavityExportCanvasPixels',[fw fh]);
for k=1:nr
    a=limits(k); levels=linspace(-a,a,cfg.contour_levels+1);
    for j=1:nc
        pos=[left+(j-1)*(side+gapx),fh-top-side-(k-1)*(side+gapy),side,side]./[fw fh fw fh];
        ax=axes('Parent',fig,'Units','normalized','Position',pos);
        draw_field(ax,Results{k},j,a,cmap,levels,cfg);
        ax.PositionConstraint='innerposition'; ax.Position=pos;
    end
    fprintf('Re=%d: shared limits for k=1/400/800 are +/- %.6g.\n',Results{k}.Re,a);
end
export_cavity_figure(fig,'fig5b',cfg);
if cfg.export_individual_fields
    for k=1:nr
        a=limits(k); levels=linspace(-a,a,cfg.contour_levels+1);
        for j=1:nc
            f=figure('Color','w','Position',[100 100 700 650],'NumberTitle','off');
            ax=axes('Parent',f); draw_field(ax,Results{k},j,a,cmap,levels,cfg);
            export_cavity_figure(f,sprintf('Re%d_snapshot_%03d',Results{k}.Re,cfg.snapshot_indices(j)),cfg);
        end
    end
end
end
function draw_field(ax,R,j,a,cmap,levels,cfg)
W=max(min(R.vorticity{j},a),-a);
contourf(ax,R.xi,R.xi,W,levels,'LineStyle','none');
axis(ax,'square'); xlim(ax,[-1 1]); ylim(ax,[-1 1]);
caxis(ax,[-a a]); colormap(ax,cmap);
set(ax,'XTick',[-1 0 1],'YTick',[-1 0 1]); xlabel(ax,'x','FontAngle','italic','FontWeight','bold'); ylabel(ax,'y','FontAngle','italic','FontWeight','bold');
style_cavity_axes(ax,cfg);
if cfg.show_panel_labels
    title(ax,['\bf\it Re\rm\bf = ' num2str(R.Re) ', \it k\rm\bf = ' num2str(R.snapshot_indices(j))], ...
        'FontName',cfg.font_name,'FontSize',cfg.panel_label_size,'FontWeight','bold','Interpreter','tex');
end
end


function style_cavity_axes(ax,cfg)
set(ax,'XTickLabelRotation',0,'FontName','Times New Roman','FontSize',cfg.font_size,'FontWeight','bold', ...
    'LineWidth',1.1,'TickDir','out','Box','on','TickLabelInterpreter','tex');
set([ax.XLabel ax.YLabel],'FontName','Times New Roman', ...
    'FontSize',cfg.label_size,'FontWeight','bold');
if cfg.grid, grid(ax,'on'); ax.GridAlpha=0.1; else, grid(ax,'off'); end
set(ax,'XMinorGrid','off','YMinorGrid','off');
end


function cmap=cavity_smooth_redblue(n)
% Continuous blue-white-red map. No central white plateau or jet end colors.
if nargin<1, n=257; end
anchors=[0.04 0.08 0.72; 0.18 0.27 0.95; 0.62 0.68 1.00; ...
         1.00 1.00 1.00; 1.00 0.65 0.65; 0.95 0.18 0.16; 0.72 0.03 0.02];
cmap=interp1(linspace(0,1,size(anchors,1)),anchors,linspace(0,1,n),'linear');
cmap=max(min(cmap,1),0);
end


function export_cavity_figure(fig,stem,cfg)
% Manuscript output only: fig5a.eps and fig5b.png.
if strcmp(stem,'fig5a')
    print(fig,fullfile(cfg.output_dir,'fig5a.eps'),'-depsc2','-painters');
elseif strcmp(stem,'fig5b')
    print(fig,fullfile(cfg.output_dir,'fig5b.png'),'-dpng',sprintf('-r%d',cfg.png_resolution));
end
end
