function [fig12a,fig12b]=plot_fig12(D)
% Figure 12: plotting only; all numerical results are supplied by main_4.
ranks=D.ranks; rhoD=D.rhoD; rhoG=D.rhoG; meanD=D.meanD; meanG=D.meanG; H=D.H;
C=taira_colors();
% ---- Figure 1: spectral radius ----
fig=figure('Name','Rank robustness - spectral radius','Color','w','Position',[100 100 820 500]);
plot(ranks,rhoD,'-o','Color',C.blue,'LineWidth',2.1,'MarkerSize',7,'MarkerFaceColor','w'); hold on;
plot(ranks,rhoG,'-s','Color',C.red,'LineWidth',2.1,'MarkerSize',7,'MarkerFaceColor','w');
yline(1,'--','Color',C.gray,'LineWidth',1.2,'HandleVisibility','off');
xlabel('\itr'); ylabel('\rho(\itA\rm)');
title('Spectral stability versus retained POD rank');
legend('DMD','GC-DMD','Location','best'); apply_taira_style(fig,'cavity19k'); ax=gca; xlim(ax,[min(ranks) max(ranks)]); xticks(ax,unique([min(ranks) 60 80 100 max(ranks)])); set(ax,'YScale','linear'); ylim(ax,[1 1.005]); yticks(ax,[1 1.0025 1.005]);

fig12a=fig;
% ---- Figure 2: long-horizon prediction error ----
fig=figure('Name','Rank robustness - long prediction','Color','w','Position',[120 120 820 500]);
semilogy(ranks,meanD,'-o','Color',C.blue,'LineWidth',2.1,'MarkerSize',7,'MarkerFaceColor','w'); hold on;
semilogy(ranks,meanG,'-s','Color',C.red,'LineWidth',2.1,'MarkerSize',7,'MarkerFaceColor','w');
xlabel('\itr'); ylabel('E_{\psi}^{int}');
set(gca,'YScale','linear');
title(sprintf('Long-horizon prediction robustness (H = %d)',H));
legend('DMD','GC-DMD','Location','best'); apply_taira_style(fig,'cavity19k'); ax=gca; xlim(ax,[min(ranks) max(ranks)]); xticks(ax,unique([min(ranks) 60 80 100 max(ranks)])); set(ax,'YScale','linear'); ylim(ax,[40 120]); yticks(ax,[40 80 120]); ax.YLabel.String='E_{\psi}^{int}'; ax.YLabel.Interpreter='tex';


fig12b=fig;
end
