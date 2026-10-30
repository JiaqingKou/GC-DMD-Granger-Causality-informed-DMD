function plot_fig3(Rgc,idx,r,operator_clim,C)
% Change combined_figure_position or the normalized positions below to tune proportions.
combined_figure_position=[30 40 1480 860];
f=figure('Name','Six lambda operators and parameter selection','NumberTitle','off', ...
    'Color','w','Position',combined_figure_position);
labels={'Unregularized','Weak regularization','Prediction optimum','GC-support optimum','Strong regularization','Over-regularized'};
marker_labels={'0','\lambda_{weak}','\lambda_{pred}','\lambda_{gc}','\lambda_{strong}','\lambda_{over}'};
map_axes=gobjects(1,6);
map_left=0.055; map_width=0.133; map_gap=0.013;
map_y=0.700; map_height=0.205; % Dedicated upper row, with separate legend band below.
for m=1:6
    map_axes(m)=axes('Parent',f,'Units','normalized', ...
        'Position',[map_left+(m-1)*(map_width+map_gap) map_y map_width map_height]);
    axm=map_axes(m); imagesc(axm,Rgc.A{idx(m)}); axis(axm,'image');
    caxis(axm,[-operator_clim operator_clim]);
    xlabel(axm,'Source j');
    if m==1, ylabel(axm,'Target i'); else, axm.YTickLabel=[]; end
end
colormap(f,taira_diverging(257));
cb=colorbar(map_axes(6)); cb.Units='normalized';
cb.Position=[0.935 map_y 0.012 map_height]; cb.Label.String='';
% Colorbar creation can shrink its peer; restore the sixth panel explicitly.
map_axes(6).Position=[map_left+5*(map_width+map_gap) map_y map_width map_height];
ax=axes('Parent',f,'Units','normalized','Position',[0.09 0.120 0.82 0.360]); hold(ax,'on');
pos=Rgc.lambda>0; lampos=Rgc.lambda(pos);
if isempty(lampos), error('Lambda grid has no positive entries.'); end
xmin=1e-7; xmax=1e-2; % Display-only decade boundaries; fitted grid unchanged.
xx=Rgc.lambda(pos);
yyaxis(ax,'left');
hRmse=semilogx(ax,xx,Rgc.val_rmse(pos),'-o','Color',C.blue,'LineWidth',2.6,'MarkerSize',5);
ylabel(ax,'RMSEval(A)');
yyaxis(ax,'right');
hF1=semilogx(ax,xx,Rgc.gc_f1(pos),'-s','Color',C.orange,'LineWidth',2.6,'MarkerSize',5);
ylabel(ax,'F1'); ylim(ax,[0 1]);
set(ax,'XScale','log'); xlim(ax,[xmin xmax]); xlabel(ax,'\lambda');

marker_colors=[C.blue; C.orange; C.yellow; 0.38 0.52 0.70; 0.20 0.55 0.28; 0.75 0.15 0.60];
marker_labels={'0','\lambda_{weak}','\lambda_{pred}','\lambda_{gc}','\lambda_{strong}','\lambda_{over}'};
% Keep marker lines separate from their labels: larger labels must not sit on curves.
for m=2:6
    k=idx(m);
    if Rgc.lambda(k)>0
        hx=xline(ax,Rgc.lambda(k),'--','Color',marker_colors(m,:),'LineWidth',1.8);
        hx.HandleVisibility='off';
    end
end
lg=legend(ax,[hRmse hF1],{'RMSEval(A)','F1'}, ...
    'Location','northoutside','Orientation','horizontal','AutoUpdate','off');
set(ax,'XScale','log'); xlim(ax,[xmin xmax]);


apply_taira_style(f);
% Reassert row-specific typography after the common style clears titles.
caption_names={'DMD',sprintf('Weak\nregularization'),sprintf('Prediction\noptimum'), ...
    sprintf('GC-support\noptimum'),sprintf('Strong\nregularization'),'Over-regularized'};
for m=1:6
    axm=map_axes(m);
    if m==1, caption='\lambda = 0';
    else, caption=sprintf('%s = %.2e',marker_labels{m},Rgc.lambda(idx(m))); end
    set(axm,'FontSize',18,'XTick',[1 5:5:r],'YTick',[1 5:5:r],'LineWidth',1.0);
    set([axm.XLabel axm.YLabel],'FontSize',20);
    title(axm,caption,'FontName','Times New Roman','FontSize',16, ...
        'FontWeight','bold','Interpreter','tex');
end
set(cb,'FontName','Times New Roman','FontSize',18,'FontWeight','bold'); cb.Label.FontSize=20;
set(lg,'FontSize',18,'Orientation','horizontal','Location','none','Units','normalized');
% Fixed legend band: neither the upper x labels nor the lower axis uses it.
lg.Position=[0.23 0.535 0.56 0.045];
ax.Position=[0.09 0.120 0.82 0.360];
set(ax,'XScale','log'); xlim(ax,[xmin xmax]);
% Include the left boundary as a labelled logarithmic tick.
xt=10.^(-7:floor(log10(xmax)));
xticks(ax,xt);
xticklabels(ax,arrayfun(@(v)sprintf('10^{%d}',round(log10(v))),xt,'UniformOutput',false));
xtickangle(ax,0);
% Normalized text coordinates work with either side of yyaxis.
% All lambda labels share one height; white backing prevents line occlusion.
label_y=0.52; % Common horizontal baseline for all five lambda labels.
for m=2:6
    v=Rgc.lambda(idx(m));
    if v<=0, continue; end
    u=(log10(v)-log10(xmin))/(log10(xmax)-log10(xmin));
    if m==3 || m==6
        u=max(0.02,u-0.012); align='right';
    else
        u=min(0.98,u+0.012); align='left';
    end
    text(ax,u,label_y,marker_labels{m},'Units','normalized', ...
        'FontName','Times New Roman','FontSize',22,'FontWeight','bold', ...
        'Interpreter','tex','Color',marker_colors(m,:), ...
        'HorizontalAlignment',align,'VerticalAlignment','middle', ...
        'BackgroundColor','w','EdgeColor','none','Margin',3, ...
        'Clipping','on','HandleVisibility','off');
end
drawnow limitrate;

end
