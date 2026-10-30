function plot_fig2(A_true,A_dmd,A_gc,A_oracle,r,operator_compare_clim)
%% Four operator figures: independently editable and exportable
ops={A_true,A_dmd,A_gc,A_oracle};
operator_names={'toy_A_true','toy_DMD','toy_GCDMD','toy_Oracle'};
for k=1:4
    fig=figure('Name',operator_names{k},'NumberTitle','off','Color','w','Position',[100 100 720 650]);
    ax=axes('Parent',fig);
    imagesc(ax,ops{k}); axis(ax,'image');
    caxis(ax,[-operator_compare_clim operator_compare_clim]);
    xlabel(ax,'Source j'); ylabel(ax,'Target i');
    colormap(fig,supp_taira_diverging(257));
    cb=colorbar(ax); cb.Label.String='';
    supp_apply_taira_style(fig);
    set(ax,'XTick',[1 5:5:r],'YTick',[1 5:5:r],'FontSize',36);
    set([ax.XLabel ax.YLabel],'FontSize',38);
    cb.FontSize=30; % No colorbar title: reserve the right margin for ticks.
    % Set the colorbar first, then restore the square peer axes.
    cb.Units='normalized'; cb.Position=[0.81 0.22 0.028 0.6092];
    ax.Units='normalized'; ax.PositionConstraint='innerposition';
    ax.Position=[0.22 0.22 0.55 0.6092];
    drawnow;
    ax.Position=[0.22 0.22 0.55 0.6092];
    % Optional: exportgraphics(fig,fullfile('figures',[operator_names{k} '.pdf']),'ContentType','vector');
end
end
