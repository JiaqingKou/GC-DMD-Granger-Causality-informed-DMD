function fig16=plot_fig16(T)
% Figure 16: plotting only. Numerical counterfactuals are computed in main_6.
names=T.support_strategy; means=T.heldout_RMSE_mean; stds=T.heldout_RMSE_std;
e_actual=T.actual_regularized_GC_RMSE(1); active_ratio=T.active_offdiag_ratio(1);
C=taira_colors(); fig=figure('Name','Equal-sparsity counterfactual test','Color','w','Position',[120 100 1440 720]); ax=axes(fig); hold(ax,'on');
b=bar(ax,1:4,means,.68,'FaceColor','flat','EdgeColor','none');
b.CData=[C.yellow;C.blue;C.green;C.lightgray];
errorbar(ax,4,means(4),stds(4),'k','LineStyle','none','LineWidth',1.3,'CapSize',10);
refline=yline(ax,e_actual,'--','Color',[.85 .325 .098],'LineWidth',2.0,'Label','Actual regularized GC-DMD','LabelHorizontalAlignment','left');
set(ax,'XTick',1:4,'XTickLabel',names); ylabel(ax,'RMSE_{test}({\itA})','Interpreter','tex');
title(ax,sprintf('Equal-sparsity counterfactual test (off-diagonal active ratio = %.1f%%)',100*active_ratio));
apply_taira_style(fig);
% Paper style for this new figure only; other figure styles are unchanged.
ax.Title.String='';
set(ax,'XGrid','off','YGrid','off','ZGrid','off','XMinorGrid','off','YMinorGrid','off','ZMinorGrid','off');
ax.XTickLabelRotation=0;
ax.TickLabelInterpreter='tex';
ax.XTickLabel={'GC-selected','Near-diagonal','Energy-local','Random'};
xlim(ax,[.5 4.5]);
% Categorical axis: all four categories are explicitly tick-labelled.
yl=ylim(ax);yticks(ax,linspace(yl(1),yl(2),5));
% Categorical axis: all four categories are explicitly tick-labelled.
yl=ylim(ax);yticks(ax,linspace(yl(1),yl(2),5));
text(ax,.125,-.18,'(refit)','Units','normalized','HorizontalAlignment','center','FontName','Times New Roman','FontSize',28,'Clipping','off'); ax.FontSize=34; ax.Position=[.19 .25 .77 .67];ax.YLabel.FontSize=36;ax.XLabel.FontSize=36;
set(ax,'FontName','Times New Roman','FontWeight','bold');
refline.FontName='Times New Roman';refline.FontSize=28;refline.FontWeight='bold';

fig16=fig;
end
