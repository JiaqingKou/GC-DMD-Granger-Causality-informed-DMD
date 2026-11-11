function [figB2a,figB2b,figB2c,figB2d] = plot_figB2(D,N,n_model,pred_horizon,dt,reconstruction_field_snapshots,prediction_field_snapshots,pod_basis_modes,error_field_snapshots,time,plot_extras,error_curve_ylim)
%PLOT_CAVITY_STORY_DIAGNOSTICS Paper-story figures for cavity GC-DMD.
% DMD=blue, GC-DMD=yellow, truth=black, POD=gray, matching project 1.
if nargin<11,plot_extras=false;end
if nargin<12,error_curve_ylim=[0 6];end
% Resolve bundled plotting dependencies before the first colormap call.
project_dir=fileparts(mfilename('fullpath'));
% Toolbox operators are provided by the unified utils directory.
C=taira_colors(); nt=D.nt;
% Match the original cavity-project vorticity rendering used in the reference figures.
% These fixed limits intentionally reproduce the original visual contrast.
field_clim = []; % shared robust scale is computed from Truth/DMD/GC-DMD in each window
error_clim = []; % errors independently scaled in each time window
cmap_vort = VorticityColormap(64);
if nargin<5 || isempty(dt), dt=0.1; end
if nargin<10 || isempty(time),time=(0:nt-1)*dt;end
time=reshape(time,1,[]);
assert(numel(time)>=nt && all(isfinite(time(1:nt))) && all(diff(time(1:nt))>0),'Invalid time vector.');
time=time(1:nt);
if nargin<6 || isempty(reconstruction_field_snapshots), reconstruction_field_snapshots=n_model; end
if nargin<7 || isempty(prediction_field_snapshots), prediction_field_snapshots=[n_model+round(pred_horizon/2), n_model+pred_horizon]; end
if nargin<8 || isempty(pod_basis_modes), pod_basis_modes=[1 3 5]; end
if nargin<9 || isempty(error_field_snapshots), error_field_snapshots=[n_model n_model+round(pred_horizon/2) nt]; end
assert(numel(error_field_snapshots)==3 && all(error_field_snapshots>=1 & error_field_snapshots<=nt),'Appendix B needs three valid snapshots.');
reconstruction_field_snapshots=unique(round(reconstruction_field_snapshots(:)'));
prediction_field_snapshots=unique(round(prediction_field_snapshots(:)'));

% Robust snapshot validation.  Snapshot numbers here are ABSOLUTE physical
% snapshot indices (not offsets from n_model).  Do not abort the whole run
% because one plotting snapshot is outside the available interval: report it
% and keep the valid requests.
rec_bad = reconstruction_field_snapshots < 1 | reconstruction_field_snapshots > n_model;
if any(rec_bad)
    warning('Ignoring reconstruction snapshots outside [1,%d]: %s', ...
        n_model, mat2str(reconstruction_field_snapshots(rec_bad)));
    reconstruction_field_snapshots = reconstruction_field_snapshots(~rec_bad);
end
pred_bad = prediction_field_snapshots <= n_model | prediction_field_snapshots > nt;
if any(pred_bad)
    warning('Ignoring prediction snapshots outside [%d,%d]: %s', ...
        n_model+1, nt, mat2str(prediction_field_snapshots(pred_bad)));
    prediction_field_snapshots = prediction_field_snapshots(~pred_bad);
end
if isempty(reconstruction_field_snapshots), reconstruction_field_snapshots=n_model; end
if isempty(prediction_field_snapshots), prediction_field_snapshots=nt; end
fprintf('Story fields: reconstruction snapshots = %s; prediction snapshots = %s (available 1:%d).\n', ...
    mat2str(reconstruction_field_snapshots), mat2str(prediction_field_snapshots), nt);

%% Recursive reconstruction/prediction from snapshot 1
ed=D.ed; eg=D.eg;
fig=figure('Name','Recursive reconstruction and prediction','Color','w','Position',[90 90 1200 650]);
figB2a=fig;
plot(time,ed,'Color',C.blue,'LineWidth',1.5); hold on; plot(time,eg,'--','Color',[0.85 0.32 0.16],'LineWidth',1.8);
ax=gca; ylim(ax,error_curve_ylim); yl=ylim(ax);
marker_y=yl(2)-.06*diff(yl);
for it=error_field_snapshots
    if it==400,continue;end % Omit only this curve marker, not field diagnostics.
    xline(ax,time(it),'--','Color',C.yellow,'LineWidth',1.5,'HandleVisibility','off');
    ha='center'; if it>=.95*nt,ha='right';end
    text(ax,time(it),marker_y,sprintf('\\it k\\rm\\bf=%d',it),'Color','k', ...
        'HorizontalAlignment',ha,'VerticalAlignment','top', ...
        'FontName','Times New Roman','FontSize',18,'FontWeight','bold', ...
        'BackgroundColor','w','Margin',1,'Clipping','on','Interpreter','tex');
end
time_tol=1e-8; % avoid roundoff adding a spurious whole time unit
left_time=floor(time(1)+time_tol); right_time=ceil(time(end)-time_tol);
if right_time<=left_time,right_time=left_time+max(dt,eps);end
xlim([left_time right_time]);
span=time(end)-time(1); mag=10^floor(log10(span/6));
choices=[1 2 5 10]*mag; tick_step=choices(find(choices>=span/6,1));
ticks=ceil(time(1)/tick_step)*tick_step:tick_step:floor(time(end)/tick_step)*tick_step;
ticks=unique([left_time ticks right_time]);
xticks(ticks);
xlabel('\itt'); ylabel('E_{\psi}(\itk)'); title('Recursive reconstruction and prediction'); legend('DMD','GC-DMD','Location','best'); apply_taira_style(fig);

%% Appendix B Figure B2(b-d): original three independent field windows
figsB2 = plot_field_set([reconstruction_field_snapshots prediction_field_snapshots], ...
    'Appendix B',D,field_clim,error_clim,cmap_vort);
assert(numel(figsB2)==3,'Appendix B requires exactly three field snapshots.');
figB2b=figsB2(1); figB2c=figsB2(2); figB2d=figsB2(3);
end

function figs = plot_field_set(indices,phase_name,D,field_clim,error_clim,cmap_vort)
% One window per time; both errors use the DMD display scale.
figs=gobjects(1,numel(indices));
fig_count=0;
for it=reshape(indices,1,[])
    xi=D.xi; yi=D.yi; FF=D.fields{it}; wT=FF{1}; wD=FF{2}; wG=FF{3}; names={'Reference','DMD','GC-DMD','DMD error','GC-DMD error'};
    physical_limits=prctile([wT(:);wD(:);wG(:)],[1 99]);
    a=prctile(abs(FF{4}(:)),99); if ~isfinite(a)||a<=0,a=1;end
    % Round outward to simple 1/2/5-based ticks; common scale from DMD.
    base=10^floor(log10(a/3)); choices=[1 2 5 10]*base;
    step=choices(find(choices>=a/3,1));
    a=ceil(a/step)*step;
    % Paper-requested display ranges at the three selected snapshots.
    switch it
        case 800, a=3; step=1;
        case 1000, a=10; step=5;
        case 1200, a=20; step=10;
    end
    shared_error_limits=[-a a];
    window_offset=40*find([800 1000 1200]==it,1);
    if isempty(window_offset),window_offset=0;end
    fig=figure('Name',sprintf('%s fields and errors k=%d',phase_name,it), ...
        'NumberTitle','off','Tag',sprintf('Cavity30k_Fields_k%d',it), ...
        'Color','w','Position',[50+window_offset 100+window_offset 1650 390]);
    fig_count=fig_count+1; figs(fig_count)=fig;
    % Uniform spacing; one shared error colorbar beside the final panel.
    xs=[.060 .237 .414 .591 .768]; aw=.137; ah=.60; y0=.18;
    for j=1:5
        ax=axes(fig,'Position',[xs(j) y0 aw ah]);
        if j<=3, limits=physical_limits; else
            limits=shared_error_limits;
        end
        draw_reference_flow(ax,xi,yi,FF{j},limits);
        % Match the 13k reference: move y label closer to its flow panel.
        ax.YLabel.Units='normalized';
        yp=ax.YLabel.Position; yp(1)=-0.105; ax.YLabel.Position=yp;
        paper_panel_label(ax,names{j});
        if j==5
            cb=colorbar(ax); cb.Position=[xs(j)+aw+.008 y0 .009 ah];
            cb.Ticks=(-round(a/step):round(a/step))*step;
            cb.TickLabels=arrayfun(@(v)sprintf('%.3g',v),cb.Ticks,'UniformOutput',false);
            cb.Label.String='\Delta\omega';
            ax.Position=[xs(j) y0 aw ah];
        end
    end
    apply_taira_style(fig);
    fprintf('Created field/error figure k=%d; shared error limits [%g,%g].\n',it,-a,a);
end
end

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
contourf(ax,x,y,Wp,linspace(limits(1),limits(2),51),'LineStyle','none');
axis(ax,'square'); xlim(ax,[-1 1]); ylim(ax,[-1 1]);
caxis(ax,limits); colormap(ax,VorticityColormap(64));
set(ax,'XTick',[-1 0 1],'YTick',[-1 0 1]); box(ax,'on'); grid(ax,'off');
xlabel(ax,'\itx'); ylabel(ax,'\ity');
end


function paper_panel_label(ax,str)
% Annotation survives the global removal of axes/figure titles.
text(ax,0.5,1.06,str,'Units','normalized','HorizontalAlignment','center', ...
 'FontName','Times New Roman','FontSize',14,'FontWeight','bold','Interpreter','tex','Clipping','off');
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
