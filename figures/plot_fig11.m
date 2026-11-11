function [fig11a,fig11b,fig11c]=plot_fig11(Md,Mg,Trec,Tpred,dt)
% Figure 11: original two Kou-truncation panels and frequency-growth map.
draw_kou_loss_19k(Trec,115,'reconstruction');fig11a=gcf;
draw_kou_loss_19k(Tpred,115,'prediction');fig11b=gcf;
C=taira_colors();wd=log(Md.lambda)/dt;wg=log(Mg.lambda)/dt;
fd=abs(imag(wd))/(2*pi);fg=abs(imag(wg))/(2*pi);gd=real(wd);gg=real(wg);
fig=figure('Name','Frequency-growth map','Color','w','Position',[80 80 620 520]); ax=axes(fig); plot(ax,fd,gd,'o','Color',C.blue,'MarkerFaceColor','w','LineWidth',1.2,'MarkerSize',6); hold(ax,'on'); plot(ax,fg,gg,'x','Color',C.red,'LineWidth',1.5,'MarkerSize',7); yline(ax,0,'k--','HandleVisibility','off'); xlim(ax,[0 1]); xticks(ax,0:.25:1); ylim(ax,[-1.2 .2]); yticks(ax,-1.2:.2:.2); xlabel(ax,'\itf'); ylabel(ax,'\sigma'); title(ax,'(b) Frequency-growth-rate map'); legend(ax,'DMD','GC-DMD','Location','best');

align_paper_panel(fig,'cavity19k'); apply_taira_style(fig,'cavity19k');
fig11c=fig;
end

function draw_kou_loss_19k(T,max_modes,interval_name)
% Original Figure 11 drawing commands; all losses supplied by main.
loss=[T.DMDLossPercent T.GCDMDLossPercent];
C=taira_colors();
fig=figure('Name',['Kou-ranked modal truncation - ' interval_name], ...
    'Color','w','Position',[100 100 900 600]);
ax=axes('Parent',fig,'Position',[.17 .17 .78 .76]);
hd=plot(ax,T.RetainedDMDModes,T.DMDLossPercent,'-o','Color',C.blue, ...
    'LineWidth',1.8,'MarkerSize',4,'MarkerFaceColor','w');hold(ax,'on');
hg=plot(ax,T.RetainedDMDModes,T.GCDMDLossPercent,'-s','Color',C.red, ...
    'LineWidth',1.8,'MarkerSize',4,'MarkerFaceColor','w');
xlabel(ax,'\itr');ylabel(ax,'E_{\psi}^{int}');
xlim(ax,[0 max_modes]);xticks(ax,unique([0:20:max_modes max_modes]));ylmax=ceil(1.08*max(loss(:))/50)*50; ylim(ax,[0 ylmax]); yticks(ax,linspace(0,ylmax,6));
legend(ax,[hd hg],{'DMD','GC-DMD'},'Location','northeast');
apply_taira_style(fig,'cavity19k');
% No MAT/CSV output: manuscript figures only.

end
