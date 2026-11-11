function [fig10a,fig10b,fig10c,fig10d]=plot_fig10(D,N,n_model,pred_horizon,dt,reconstruction_field_snapshots,prediction_field_snapshots,pod_basis_modes,error_field_snapshots,time)
% Figure 10: original recursive error and three field-comparison panels.
C=taira_colors(); nt=D.nt;
field_clim=[]; error_clim=[]; cmap_vort=VorticityColormap(64);
time=reshape(time,1,[]); time=time(1:nt);
% Respect the exact snapshot indices configured in main_4 (no automatic fourth slice).
error_field_snapshots=round(error_field_snapshots(:)');
assert(numel(error_field_snapshots)==3 && all(isfinite(error_field_snapshots)) && ...
    all(error_field_snapshots>=1 & error_field_snapshots<=nt) && ...
    numel(unique(error_field_snapshots))==3, ...
    'Figure 10 requires three distinct error-field snapshots within 1:%d.',nt);

ed=D.ed; eg=D.eg;
fig=figure('Name','Recursive reconstruction and prediction','Color','w','Position',[90 90 1200 650]);
plot(time,ed,'Color',C.blue,'LineWidth',1.5); hold on; plot(time,eg,'--','Color',[0.85 0.32 0.16],'LineWidth',1.8);
ax=gca; yl=ylim(ax); ylim(ax,[yl(1) yl(2)+.10*diff(yl)]);
for it=error_field_snapshots
    if it==400,continue;end
    xline(ax,time(it),'--','Color',C.yellow,'LineWidth',1.5,'HandleVisibility','off');
end

xlim([floor(time(1)) time(end)]);
span=time(end)-time(1); mag=10^floor(log10(span/6));
choices=[1 2 5 10]*mag; tick_step=choices(find(choices>=span/6,1));
ticks=ceil(time(1)/tick_step)*tick_step:tick_step:floor(time(end)/tick_step)*tick_step;
xticks(unique([floor(time(1)) ticks time(end)])); yl=ylim(ax); top=ceil(yl(2)); ylim(ax,[0 top]); yticks(ax,linspace(0,top,6));
xlabel('\itt'); ylabel('E_{\psi}(\itk\rm)'); title('Recursive reconstruction and prediction'); legend('DMD','GC-DMD','Location','best'); apply_taira_style(fig,'cavity19k');
% Draw labels last, after axis limits and style; preserve original alignment.
ax=gca; lim=ylim(ax); label_y=lim(2)-0.12*diff(lim);
for it=error_field_snapshots
    if it==400,continue;end
    ha='center'; if it>=.95*nt,ha='right';end
    label=['{\it k}=' num2str(it)];
    text(ax,time(it),label_y,label,'Interpreter','tex', ...
        'FontName','Times New Roman','FontSize',18,'FontWeight','bold', ...
        'Color','k','HorizontalAlignment',ha,'VerticalAlignment','middle', ...
        'BackgroundColor','white','Margin',2,'Clipping','off');
end



fig10a=fig;
field_figs=[plot_field_set(reconstruction_field_snapshots,'Reconstruction',D,field_clim,error_clim,cmap_vort), ...
    plot_field_set(prediction_field_snapshots,'Prediction',D,field_clim,error_clim,cmap_vort)];
assert(numel(field_figs)==3,'Figure 10 requires three snapshot windows.');
fig10b=field_figs(1);fig10c=field_figs(2);fig10d=field_figs(3);
end

function figs=plot_field_set(indices,phase_name,D,field_clim,error_clim,cmap_vort)
% One window per time; both errors use the DMD display scale.
figs=gobjects(1,numel(indices));
for it=indices
    xi=D.xi; yi=D.yi; FF=D.fields{it}; wT=FF{1}; wD=FF{2}; wG=FF{3}; names={'Reference','DMD','GC-DMD','DMD error','GC-DMD error'};
    physical_limits=prctile([wT(:);wD(:);wG(:)],[1 99]);
    a=prctile(abs(FF{4}(:)),99); if ~isfinite(a)||a<=0,a=1;end
    % Round outward to simple 1/2/5-based ticks; common scale from DMD.
    base=10^floor(log10(a/3)); choices=[1 2 5 10]*base;
    step=choices(find(choices>=a/3,1));
    a=ceil(a/step)*step;
    % Paper-requested display ranges at the three selected snapshots.
    switch it
        case 800, a=0.4; step=0.2;
        case 1000, a=1.6; step=0.8;
        case 1200, a=4; step=2;
    end
    shared_error_limits=[-a a];
    fig=figure('Name',sprintf('%s fields and errors t=%d',phase_name,it),'Color','w','Position',[50 100 1650 390]);
    % Uniform spacing; one shared error colorbar beside the final panel.
    xs=[.060 .237 .414 .591 .768]; aw=.137; ah=.60; y0=.18;
    for j=1:5
        ax=axes(fig,'Position',[xs(j) y0 aw ah]);
        if j<=3, limits=physical_limits; else
            limits=shared_error_limits;
        end
        draw_reference_flow(ax,xi,yi,FF{j},limits);
        % Move y labels closer to their respective flow panels.
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
    apply_taira_style(fig,'cavity19k');
    figs(find(indices==it,1,'first'))=fig;
end
end

function paper_panel_label(ax,str)
% Annotation survives the global removal of axes/figure titles.
text(ax,0.5,1.06,str,'Units','normalized','HorizontalAlignment','center', ...
 'FontName','Times New Roman','FontSize',14,'FontWeight','bold','Interpreter','tex','Clipping','off');
end
