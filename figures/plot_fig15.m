function [fig15a,fig15b]=plot_fig15(R,P)
% Figure 15(a,b): two original plotting windows, unchanged.
plot_lambda_generalization(R); fig15a=gcf;
plot_gc_pareto(P); fig15b=gcf;
end

function plot_lambda_generalization(R)
% Log scale explicitly set before hold; mark the actually selected lambda.
C=taira_colors(); lam=R.lambda_grid(:); err=R.cv_gc(:);
valid=lam>0 & isfinite(lam) & isfinite(err); L=lam(valid); E=err(valid);
assert(~isempty(L),'No valid positive lambda values for plotting.');
selected=R.lambda_gc;
fig=figure('Name','Regularization selection','Color','w','Position',[100 100 1360 650]);
ax=axes(fig);set(ax,'XScale','log');hold(ax,'on');
h=semilogx(ax,L,E,'-s','Color',C.red,'MarkerFaceColor','w','MarkerSize',4,'LineWidth',2.0);
hd=yline(ax,R.cv_dmd,'--','Color',C.gray,'LineWidth',1.8);
if selected>0 && isfinite(selected)
 xline(ax,selected,'--','Color',C.blue,'LineWidth',1.4,'HandleVisibility','off');
 [~,ib]=min(abs(lam-selected));
 if isfinite(err(ib)),plot(ax,selected,err(ib),'s','Color',C.red,'MarkerFaceColor',C.red,'MarkerSize',8,'HandleVisibility','off');end
 center=selected;
else
 [~,ib]=min(E);center=L(ib);
 text(ax,.03,.07,'\lambda_{pred} = 0','Units','normalized','Interpreter','tex','Color',C.blue);
end
lo=1e-8;hi=1e-3;
xlim(ax,[lo hi]);
xlabel(ax,'{\it \lambda}','Interpreter','tex');ylabel(ax,'RMSE_{val}({\itA})','Interpreter','tex');
legend(ax,[h hd],{'GC-DMD','DMD (\lambda=0)'},'Location','best');
apply_taira_style(fig);ax.Title.String='';grid(ax,'off');
set(ax,'XGrid','off','YGrid','off','XMinorGrid','off','YMinorGrid','off','XMinorTick','on');
% Show every decade plus 2:9 minor ticks; minor grid remains off.
xl=xlim(ax);decades=floor(log10(xl(1))):ceil(log10(xl(2)));
major=10.^decades;xticks(ax,major(major>=xl(1) & major<=xl(2)));
% Include lower/upper endpoints and use readable evenly spaced y ticks.
yl=ylim(ax);yticks(ax,linspace(yl(1),yl(2),5));
% Include lower/upper endpoints and use readable evenly spaced y ticks.
yl=ylim(ax);yticks(ax,linspace(yl(1),yl(2),5));
if isprop(ax.XAxis,'MinorTickValues')
 minor=reshape((2:9)'*(10.^decades),1,[]);
 ax.XAxis.MinorTickValues=sort(minor(minor>xl(1) & minor<xl(2)));
end
if selected>0 && isfinite(selected)
 yl=ylim(ax);tx=selected*10^(.015*diff(log10(xl)));
 text(ax,tx,yl(1)+.72*diff(yl),'\lambda_{pred}','Rotation',90, ...
  'Interpreter','tex','FontName','Times New Roman','FontSize',30,'FontWeight','bold', ...
  'Color',C.blue,'HorizontalAlignment','center','VerticalAlignment','top', ...
  'BackgroundColor','w','Margin',1,'Clipping','on');
end
end

function plot_gc_pareto(P)
C=taira_colors(); fig=figure('Color','w','Position',[100 100 1250 650]); ax=axes(fig); hold(ax,'on');
plot(ax,P.active_gc,P.pred_gc,'-s','Color',C.red,'MarkerFaceColor','w','LineWidth',2.0); plot(ax,P.active_dmd,P.pred_dmd,'o','Color',C.blue,'MarkerFaceColor',C.blue,'MarkerSize',8);
if ~isempty(P.selected_idx), k=P.selected_idx; plot(ax,P.active_gc(k),P.pred_gc(k),'s','Color',C.red,'MarkerFaceColor',C.red,'MarkerSize',9,'HandleVisibility','off');end
xlabel(ax,'NNZ({\itA})/{\it r}^2'); ylabel(ax,'RMSE_{test}({\itA})'); title(ax,'Sparsity--generalization Pareto curve'); legend(ax,'GC-DMD','DMD','Location','best'); apply_taira_style(fig);ax.YLabel.FontSize=30;ax.XLabel.FontSize=32;ax.FontSize=30;ax.Position=[.23 .24 .72 .65];
xl=xlim(ax);yl=ylim(ax);xticks(ax,linspace(xl(1),xl(2),6));yticks(ax,linspace(yl(1),yl(2),5));
xl=xlim(ax);yl=ylim(ax);xticks(ax,linspace(xl(1),xl(2),6));yticks(ax,linspace(yl(1),yl(2),5));
end
