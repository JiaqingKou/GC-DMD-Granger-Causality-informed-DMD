function fig=plot_fig14(G,Ad,Ag,R)
% Reference cavity geometry; axes belong directly to the figure.
% Panel letters are a second xlabel line, never outside a clipping panel.
R=R(R.rank>=100,:);blue=[0 .447 .741];red=[.85 .325 .098];
fig=figure('Name','GC operators and rank sparsity','Color','w', ...
 'Units','pixels','Position',[50 40 1100 1100],'WindowStyle','normal');
movegui(fig,'center');
pp={[.10 .60 .30 .30],[.59 .60 .30 .30],[.10 .12 .30 .30],[.59 .12 .30 .30]};
axs=gobjects(4,1);cbs=gobjects(3,1);
off=~eye(size(Ad));v=[abs(Ad(off));abs(Ag(off))];a=prctile(v,99.5);if ~isfinite(a)||a<=0,a=1;end
v=G(off);v=v(isfinite(v)&v>0);if isempty(v),gmax=1;else,gmax=prctile(v,99);end
mats={G,Ad,Ag};
for k=1:4
 ax=axes('Parent',fig,'Units','normalized','Position',pp{k});axs(k)=ax;
 set(ax,'FontName','Times New Roman','FontSize',22,'FontWeight','bold','Box','on', ...
  'XGrid','off','YGrid','off','XMinorGrid','off','YMinorGrid','off','TickDir','out');
 if k<=3
  imagesc(ax,mats{k});ax.YDir='reverse';r=size(G,1);xlim(ax,[.5 r+.5]);ylim(ax,[.5 r+.5]);
  ticks=unique(round(linspace(1,r,3)));xticks(ax,ticks);yticks(ax,ticks);
  if k==1
   u=linspace(0,1,256)';colormap(ax,(1-u)*[1 1 1]+u*[.72 .12 .12]);caxis(ax,[0 gmax]);
   xlabel(ax,'Cause {\it j}');ylabel(ax,'Effect {\it i}');tag='GC';
  else
   colormap(ax,taira_diverging(257));caxis(ax,[-a a]);xlabel(ax,'Source {\it j}');ylabel(ax,'Target {\it i}');tag='A_{ij}';
  end
  cb=colorbar(ax);cbs(k)=cb;cb.Units='normalized';
  cb.Ticks=linspace(ax.CLim(1),ax.CLim(2),5);
  cb.TickLabels=arrayfun(@(x)sprintf('%.2f',x),cb.Ticks,'UniformOutput',false);
  cb.FontName='Times New Roman';cb.FontWeight='bold';
  cb.Title.String='';cb.Title.FontName='Times New Roman';cb.Title.Interpreter='tex';
 else
  plot(ax,R.rank,R.active_DMD,'-o','Color',blue,'LineWidth',1.5);hold(ax,'on');
  plot(ax,R.rank,R.active_GC,'-s','Color',red,'LineWidth',1.5);
  xlim(ax,[min(R.rank) max(R.rank)]);xticks(ax,unique([min(R.rank) 200 max(R.rank)]));ylim(ax,[0 1]);yticks(ax,[0 .2 .4 .6 .8 1]);
  xlabel(ax,'{\it r}');ylabel(ax,'NNZ({\it A})/{\it r}^2');
  legend(ax,'DMD','GC-DMD','FontName','Times New Roman','Box','off','Location','best');
 end
 ax.PositionConstraint='innerposition';ax.DataAspectRatioMode='auto';
 ax.PlotBoxAspectRatio=[1 1 1];ax.PlotBoxAspectRatioMode='manual';
 % Match the reference project: reserve a second xlabel line for (a)--(d).
 axis_label=ax.XLabel.String;
 ax.XLabel.String={axis_label,sprintf('(%c)','a'+k-1)};
 set([ax.XLabel ax.YLabel],'FontName','Times New Roman', ...
  'FontSize',24,'FontWeight','bold','Interpreter','tex');
end
for k=1:4
 ax=axs(k);ax.PositionConstraint='innerposition';ax.Position=pp{k};
 % Apply AFTER imagesc/plot, which can reset axes typography.
 set(ax,'FontName','Times New Roman','FontUnits','points', ...
  'FontSize',26,'FontWeight','bold','TickLabelInterpreter','tex', ...
  'XTickLabelRotation',0,'YTickLabelRotation',0);
 % Preserve the existing axis-title sizes and layout.
 set([ax.XLabel ax.YLabel],'FontName','Times New Roman', ...
  'FontSize',24,'FontWeight','bold');
 pbaspect(ax,[1 1 1]);
 if k<=3
  cb=cbs(k);cb.Units='normalized';
  cb.Position=[pp{k}(1)+pp{k}(3)+.0075 pp{k}(2) .012 pp{k}(4)];
  cb.FontSize=20;cb.Title.FontSize=22;
  ax.Position=pp{k};
 else
  lg=findall(fig,'Type','legend');
  set(lg,'FontSize',21,'FontName','Times New Roman','FontWeight','bold');
 end
end
% Preserve geometry and point sizes during optional PNG export. No EPS export.
set(fig,'PaperPositionMode','auto');
end
