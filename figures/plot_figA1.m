function figs=plot_figA1(R,Atrue,operator_color_max)
% Figure A1(a-d): original standalone annotated operator plots.
if nargin<3,operator_color_max=.20;end
ops={Atrue,R.DMD,R.A{R.best(1),1},R.A{R.best(2),2}};
names={'True A','DMD','L1-DMD','GC-DMD'};
lim=max(cellfun(@(A)max(abs(A(:))),ops));
figs=gobjects(1,4);
for j=1:4
 f=figure('Name',names{j},'Color','w','Position',[80+25*j 80 850 850]);
 ax=axes('Parent',f,'Units','normalized','Position',[.20 .18 .73 .73]);
 imagesc(ax,ops{j});axis(ax,'image');caxis(ax,[-lim lim]);
 set(ax,'XTick',1:5,'YTick',1:5);xlabel(ax,'Source {\it j}','Interpreter','tex');ylabel(ax,'Target {\it i}','Interpreter','tex');
 title(ax,'');
 for ii=1:5
  for jj=1:5
   v=ops{j}(ii,jj);cc='k';if abs(v)>.55*lim,cc='w';end
   if abs(v)<1e-10,label='0';else,label=sprintf('%.2f',v);end
   text(ax,jj,ii,label,'HorizontalAlignment','center','FontName','Times New Roman', ...
    'FontSize',24,'FontWeight','bold','Color',cc);
  end
 end
 colormap(f,taira_diverging(257));paper_style(f);
 figs(j)=f;
end
end
function paper_style(f)
for ax=findall(f,'Type','axes')'
 set(ax,'FontName','Times New Roman','FontSize',36,'FontWeight','bold','XGrid','off','YGrid','off', ...
 'XMinorGrid','off','YMinorGrid','off','XTickLabelRotation',0,'Box','on');
 ax.XLabel.FontSize=38;ax.YLabel.FontSize=38;ax.Title.FontSize=36;ax.Title.String='';
 % Preserve each label's existing Interpreter (LaTeX or TeX).
 if isprop(ax,'Toolbar'),ax.Toolbar.Visible='off';end
end
set(findall(f,'Type','legend'),'FontName','Times New Roman','FontSize',30,'FontWeight','bold','Box','off');
for cb=findall(f,'Type','colorbar')',cb.FontName='Times New Roman';cb.FontSize=30;cb.Label.FontSize=32;end
end
