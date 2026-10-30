function plot_fig4(R,A_true,idx_show,support_tol)
%PLOT_LAMBDA_COEFFICIENT_SURVIVAL Ranked coefficient distributions + recovery error.
% Every operator is sorted independently by its own off-diagonal |A_ij|.
if nargin<4, support_tol=1e-3; end
C=supp_taira_colors(); r=size(A_true,1); mask=~eye(r);
atrue=sort(abs(A_true(mask)),'descend');
cols=[0 0 0; C.blue; C.orange; C.yellow; 0.38 0.52 0.70; 0.20 0.55 0.28; 0.75 0.15 0.60];

fig=figure('Visible','off','Name','Ranked coefficient distributions and operator error', ...
    'Color','w','Position',[70 60 1400 800]);
stage_names={'DMD','Weak regularization','Prediction optimum','GC-support optimum','Strong regularization','Over-regularized'};

%% (a) independently ranked coefficient distributions
ax=axes(fig); hold(ax,'on');
semilogy(ax,1:numel(atrue),max(atrue,1e-8),'k-','LineWidth',2.6);
leg={'True operator'};
errs=zeros(numel(idx_show),1);
for m=1:numel(idx_show)
    k=idx_show(m);
    avec=sort(abs(R.A{k}(mask)),'descend');
    cc=cols(min(m+1,size(cols,1)),:);
    semilogy(ax,1:numel(avec),max(avec,1e-8),'-','Color',cc,'LineWidth',2.6);
    errs(m)=norm(R.A{k}-A_true,'fro')/max(norm(A_true,'fro'),eps);
    leg{end+1}=stage_names{m}; %#ok<AGROW>
end
yline(ax,support_tol,'--','Color',[0.48 0.48 0.48],'LineWidth',1.3);
leg{end+1}='Support threshold';
set(ax,'YScale','log');
xlim(ax,[0 400]); xticks(ax,0:100:400);
ymax=max(10,max(atrue)*100); % Reserve top logarithmic space for the complete legend.
ymin=max(1e-5,support_tol/20);
ylim(ax,[1e-4 1e1]); yticks(ax,10.^(-4:1));
xlabel(ax,'p');
ylabel(ax,'|A_{ij}^{true}|');
title(ax,'(a) Ranked coefficient distributions');
lg=legend(ax,leg,'Location','northeast','NumColumns',1,'AutoUpdate','off');
supp_apply_taira_style(fig); set(ax,'YScale','log','YMinorTick','on');
set(ax,'FontSize',28); set([ax.XLabel ax.YLabel],'FontSize',30);
ax.YLabel.FontSize=28;
ax.Units='normalized'; ax.Position=[0.20 0.16 0.75 0.76];
set(lg,'FontSize',25,'NumColumns',2,'Location','northeast','Units','normalized');
lg.Position=[0.23 0.75 0.73 0.17];
set(lg,'FontName','Times New Roman','FontWeight','bold','Interpreter','tex');
lg.ItemTokenSize=[24 18]; set(fig,'Visible','on'); drawnow limitrate;

end
