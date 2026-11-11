function [fig7a,fig7b] = plot_fig7(D,n_model,pred_horizon,dt,reconstruction_field_snapshots,prediction_field_snapshots,pod_basis_modes,error_field_snapshots,time)
%PLOT_CAVITY_STORY_DIAGNOSTICS Paper-story figures for cavity GC-DMD.
% DMD=blue, GC-DMD=red, truth=black, POD=gray, matching project 1.
C=taira_colors(); nt=D.nt;
% Match the original cavity-project vorticity rendering used in the reference figures.
% These fixed limits intentionally reproduce the original visual contrast.
field_clim = []; % shared robust scale is computed from Truth/DMD/GC-DMD in each window
error_clim = [-0.05 0.05]; % fixed paper comparison scale for all eight error slices
cmap_vort = VorticityColormap(64);
if nargin<4 || isempty(dt), dt=0.1; end
if nargin<9 || isempty(time), time=(0:nt-1)*dt; end
time=double(time(:)'); assert(numel(time)>=nt,'Time vector is shorter than the data horizon.');
time=time(1:nt);
if nargin<5 || isempty(reconstruction_field_snapshots), reconstruction_field_snapshots=n_model; end
if nargin<6 || isempty(prediction_field_snapshots), prediction_field_snapshots=[n_model+round(pred_horizon/2), n_model+pred_horizon]; end
if nargin<7 || isempty(pod_basis_modes), pod_basis_modes=[1 3 5]; end
if nargin<8 || isempty(error_field_snapshots), error_field_snapshots=[round(n_model/2) n_model n_model+round(pred_horizon/2) nt]; end
assert(numel(error_field_snapshots)==4 && all(error_field_snapshots>=1 & error_field_snapshots<=nt),'Figure 7 needs four valid precomputed snapshots.');
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

%% Precomputed recursive prediction errors (from main)
ed=D.ed; eg=D.eg;
fig=figure('Name','Recursive reconstruction and prediction','Color','w','Position',[90 90 1300 780]);
ax=axes('Parent',fig,'Position',[0.16 0.18 0.79 0.73]);
hd=plot(ax,time,ed,'Color',C.blue,'LineWidth',2.0); hold(ax,'on');
hg=plot(ax,time,eg,'--','Color',C.red,'LineWidth',2.2);
ymax=max([ed(:);eg(:)]); if ~isfinite(ymax)||ymax<=0, ymax=1; end
ylim(ax,[0 0.05]); yticks(ax,0:0.01:0.05); % Explicit endpoints and uniform 0.01 spacing.
tleft=floor(time(1)); tright=ceil(time(end));
if tright<=tleft,tright=tleft+dt;end
xlim(ax,[tleft tright]);
span=tright-tleft; base=10^floor(log10(span/5));
choices=base*[1 2 5 10]; step=choices(find(choices>=span/6,1));
ticks=ceil(tleft/step)*step:step:floor(tright/step)*step;
ticks=unique([tleft ticks tright]);
% Avoid near-duplicate endpoint tick labels.
if numel(ticks)>2 && ticks(2)-ticks(1)<0.35*step,ticks(2)=[];end
if numel(ticks)>2 && ticks(end)-ticks(end-1)<0.35*step,ticks(end-1)=[];end
xticks(ax,ticks);
for im=1:numel(error_field_snapshots)
    it=error_field_snapshots(im); tm=time(it);
    xline(ax,tm,'--','Color',C.red,'LineWidth',1.6,'HandleVisibility','off');
    ha='center'; if im==1,ha='right';elseif im==numel(error_field_snapshots),ha='right';end
    if im==1, tm=tm-0.008*(tright-tleft); end
    text(ax,tm,ymax*1.12,sprintf('\\it k\\rm\\bf = %d',it),'Color','k', ...
        'HorizontalAlignment',ha,'VerticalAlignment','middle', ...
        'FontName','Times New Roman','FontSize',18,'FontWeight','bold', ...
        'BackgroundColor','w','Margin',1,'Clipping','on','Interpreter','tex');
end
xlabel(ax,'\itt','Interpreter','tex'); ylabel(ax,'\itE_{\psi}(k)','Interpreter','tex');
legend(ax,[hd hg],{'DMD','GC-DMD'},'Location','southeast'); apply_taira_style(fig,'cavity13k');

fig7a=fig;

% Eight signed vorticity error slices, with a common robust scale and levels.
fprintf('Preparing 2x4 error fields at k=%s...\n',mat2str(error_field_snapshots));
EE=cell(2,4); av=[];
xi=D.xi; yi=D.yi;
for col=1:4
    it=error_field_snapshots(col);
    FF=D.fields{it};
    assert(~isempty(FF),'Missing precomputed vorticity at k=%d',it);
    EE{1,col}=FF{4}; EE{2,col}=FF{5};
    av=[av;abs(EE{1,col}(:));abs(EE{2,col}(:))]; %#ok<AGROW>
end
av=av(isfinite(av));
if isempty(error_clim)
    ae=prctile(av,97.0); if isempty(ae)||~isfinite(ae)||ae<=0, ae=1; end
    ecl=[-ae ae];
else
    ecl=error_clim; ae=max(abs(ecl));
end
levelsE=linspace(-ae,ae,81);
fig=figure('Name','Signed vorticity error slices','Color','w','Position',[80 80 1440 700]);
tl=tiledlayout(fig,2,4,'TileSpacing','compact','Padding','compact');
tl.Units='normalized'; tl.OuterPosition=[0.085 0.02 0.91 0.96];
method_names={'DMD','GC-DMD'};
add_vertical_method_labels(fig,method_names);
for row=1:2
    for col=1:4
        ax=nexttile(tl,(row-1)*4+col); W=EE{row,col};
        % Original filled-contour rendering: 80 automatically selected levels.
        contourf(ax,xi,yi,W,80,'LineStyle','none'); axis(ax,'square'); xlim(ax,[-1 1]); ylim(ax,[-1 1]); box(ax,'on');
        caxis(ax,ecl); colormap(ax,cmap_vort); xlabel(ax,'\itx','Interpreter','tex'); ylabel(ax,'\ity','Interpreter','tex');
        if row==1, text(ax,0.5,1.10,sprintf('\\it k\\rm\\bf = %d',error_field_snapshots(col)),'Units','normalized','HorizontalAlignment','center','FontName','Times New Roman','FontSize',18,'FontWeight','bold','Clipping','off','Interpreter','tex'); end
    end
end
cb=colorbar(ax); cb.Layout.Tile='east'; cb.Ticks=[-0.05 0 0.05]; cb.Label.String='\Delta\omega'; apply_taira_style(fig,'cavity13k');

fprintf('2x4 error fields completed.\n');
fig7b=fig;
end % plot_fig7


function add_vertical_method_labels(fig,names)
% Dedicated in-canvas strip: vertical labels do not extend beyond the figure.
ax=axes('Parent',fig,'Units','normalized','Position',[0.055 0.02 0.035 0.96], ...
    'Visible','off','XLim',[0 1],'YLim',[0 1],'Tag','PaperRowLabels','HandleVisibility','off');
method_font=24;
if strcmp(fig.Name,'Signed vorticity error slices'),method_font=24;end
for row=1:numel(names)
    text(ax,0.5,1-(row-0.5)/numel(names),names{row}, ...
        'Rotation',90,'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'FontName','Times New Roman','FontSize',method_font,'FontWeight','bold','Clipping','on');
end
end

