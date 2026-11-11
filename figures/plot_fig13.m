function fig=plot_fig13(X,basis,nx,ny,idx)
if nargin<6 || isempty(idx),idx=[1 320 640 800];end
idx=idx(idx>=1 & idx<=size(X,2));
fields=X(:,idx); % caller supplies mean-subtracted pressure snapshots
% Dedicated plot: no generic style function can restore removed x labels.
a=prctile(abs(fields(:)),99.5);if ~isfinite(a)||a<=0,a=max(abs(fields(:)));end
n=numel(idx);pw=760;ph=pw/5;left=95;bottom=100;gap=28;top=25;
height=bottom+n*ph+(n-1)*gap+top;width=left+pw+8+16+155;
fig=figure('Name','Pressure snapshots','Color','w','Units','pixels', ...
 'Position',[50 50 width height],'WindowStyle','normal');
axs=gobjects(n,1);
for k=1:n
 ax=axes('Parent',fig,'Units','pixels','Position',[left bottom+(n-k)*(ph+gap) pw ph]);axs(k)=ax;
 imagesc(ax,linspace(0,20,ny),linspace(0,4,nx),reshape(real(fields(:,k)),nx,ny));
 set(ax,'YDir','normal','FontName','Times New Roman','FontSize',24,'FontWeight','bold', ...
  'Box','on','XGrid','off','YGrid','off','XMinorGrid','off','YMinorGrid','off','TickDir','out');
 xlim(ax,[0 20]);ylim(ax,[0 4]);yticks(ax,[0 2 4]);
 colormap(ax,taira_diverging(257));caxis(ax,[-a a]);
 ylabel(ax,'{\it r}/{\it D}','FontSize',26,'Interpreter','tex');
 if k<n
  ax.XTick=[0 5 10 15 20];ax.XTickLabel={};ax.XLabel.String='';ax.XColor='k';ax.Box='on';
 else
  xticks(ax,[0 5 10 15 20]);xlabel(ax,'{\it x}/{\it D}','FontSize',26,'Interpreter','tex');
 end
 text(ax,.02,.90,['{\it k} = ' num2str(idx(k))],'Units','normalized', ...
  'FontName','Times New Roman','FontSize',22,'FontWeight','bold','Color','k', ...
  'VerticalAlignment','top','BackgroundColor','w','Margin',1);
end
% Full-height colorbar, located eight pixels from the image edge.
cb=colorbar(axs(end));set(cb,'Units','pixels','Position',[left+pw+8 bottom 16 n*ph+(n-1)*gap], ...
 'FontName','Times New Roman','FontSize',22,'FontWeight','bold');
cb.Label.String=char([112 771]);cb.Label.Interpreter='none';cb.Label.FontWeight='bold';cb.Label.FontAngle='italic';cb.Label.FontName='Times New Roman';cb.Label.FontSize=24;
for k=1:n,axs(k).Position=[left bottom+(n-k)*(ph+gap) pw ph];end
set(fig,'PaperPositionMode','auto');
end
