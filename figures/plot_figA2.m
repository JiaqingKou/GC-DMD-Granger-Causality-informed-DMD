function f=plot_figA2(R)
% Figure A2: original lambda validation panel; no numerical changes.
blue=[0 .447 .741];red=[.85 .325 .098];
f=figure('Name','Lambda validation','Color','w','Position',[90 90 1250 850]);
ax=axes(f,'Position',[.18 .18 .78 .74]);hold(ax,'on');pos=R.lambda>0;
h1=semilogx(ax,R.lambda(pos),R.val(pos,1),'-o','Color',blue,'LineWidth',2,'MarkerSize',5);
h2=semilogx(ax,R.lambda(pos),R.val(pos,2),'-s','Color',red,'LineWidth',2,'MarkerSize',5);
h0=yline(ax,R.val(1,1),'--','Color',[.4 .4 .4]);set(ax,'XScale','log');
xlabel(ax,'\lambda','Interpreter','tex');
% Use the same MATLAB TeX markup as the working matrix labels.
% No LaTeX commands, and no literal space before italic A.
ylabel(ax,'RMSE_{val}({\itA})','Interpreter','tex', ...
 'FontName','Times New Roman','FontWeight','bold');
% Explicit endpoints and decade ticks; the first small-lambda point may be clipped.
xlim(ax,[1e-7 1e-2]);
xticks(ax,10.^(-7:-2));
% Vertical legend in the free upper-left corner.
legend(ax,[h1 h2 h0],{'L1-DMD','GC-DMD','DMD'}, ...
 'Location','northwest','Orientation','vertical','AutoUpdate','off');
% Tight, regular ticks with values at both axis endpoints.
yLo=0.095;yHi=0.120;yStep=0.005;
ylim(ax,[yLo yHi]);
yticks(ax,yLo:yStep:yHi);
box(ax,'on');
set(ax,'Layer','top','TickDir','in');
yl=ylim(ax);
for j=1:2
 k=R.best(j);col=blue;if j==2,col=red;end
 if R.lambda(k)>0
  xline(ax,R.lambda(k),'--','Color',col,'LineWidth',1.6,'HandleVisibility','off');
  plot(ax,R.lambda(k),R.val(k,j),'o','Color',col,'MarkerSize',11,'LineWidth',2,'HandleVisibility','off');
  % Place each label immediately to the RIGHT of its line, inside the axes.
  % Both labels share the red label height; blue left, red right.
  tx=R.lambda(k)*10^(0.018*5);ty=yl(1)+0.12*diff(yl);
  ha='left';if j==1,tx=R.lambda(k)/10^(0.018*5);ha='right';end
  text(ax,tx,ty,'\lambda_{pred}','Interpreter','tex','Color',col, ...
   'FontName','Times New Roman','FontSize',30,'FontWeight','bold', ...
   'HorizontalAlignment',ha,'VerticalAlignment','middle', ...
   'BackgroundColor','w','Margin',2,'Clipping','on','HandleVisibility','off');
 end
end
paper_style(f);
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
