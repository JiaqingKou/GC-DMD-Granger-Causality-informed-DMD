function plot_fig1(A_true,GC_raw,G,A_graph,r,C,true_A_tau,gc_tail_level,output_dir)
f=figure('Color','w','Position',[70 70 1400 700]); ax=axes(f); hold(ax,'on');
off=~eye(r); true_vals=abs(A_true(off)); gc_vals=GC_raw(off);
[as,ord]=sort(true_vals,'descend'); gc_paired=gc_vals(ord)/G.scale;
x=1:numel(as);
yyaxis(ax,'left'); hA=plot(ax,x,as,'-','Color',C.blue,'LineWidth',2.8);
ylabel(ax,'|A_{ij}^{true}|'); ylim(ax,[0 .5]); yticks(ax,0:.1:.5);
yyaxis(ax,'right'); hG=plot(ax,x,gc_paired,'-o','Color',C.orange,'LineWidth',2.0,'MarkerSize',3.6);
ylabel(ax,'Gtilde_{ij}'); ylim(ax,[0 1]); yticks(ax,0:.2:1);
xlabel(ax,'p');
apply_taira_style(f);
set(ax,'FontName','Times New Roman','FontSize',28,'FontWeight','bold','LineWidth',1.3);
set([ax.XLabel ax.YLabel],'FontName','Times New Roman','FontSize',30,'FontWeight','bold');
ax.Units='normalized'; ax.Position=[0.14 0.19 0.70 0.75];
xlim(ax,[0 400]); xticks(ax,0:50:400);

% Separate numeric export allows assessment of GC tail against true-A rank.
paired_table=table(x',as,gc_paired,as>true_A_tau,gc_paired>gc_tail_level, ...
 'VariableNames',{'Rank_by_true_A','Abs_true_A','Paired_GC_normalized','True_A_retained','GC_retained'});
writetable(paired_table,fullfile(output_dir,'Fig1_paired_rank_data.csv'));
% FIG 1C-1: true-A support, plotted as a standalone panel for manuscript insertion.
f=figure('Color','w','Position',[100 100 720 650]); ax=axes(f);
imagesc(ax,A_graph); axis(ax,'image'); caxis(ax,[0 1]);
xlabel(ax,'Source j'); ylabel(ax,'Target i');
set(ax,'XTick',[1 5:5:r],'YTick',[1 5:5:r]);
colormap(f,[1 1 1; 0 0 0]); % 0: white (absent), 1: black (present). apply_taira_style(f);
set(ax,'FontName','Times New Roman','FontSize',26,'FontWeight','bold','LineWidth',1.3);
set([ax.XLabel ax.YLabel],'FontName','Times New Roman','FontSize',28,'FontWeight','bold');

% FIG 1C-2: BH-FDR significant conditional-GC support, standalone panel.
f=figure('Color','w','Position',[100 100 720 650]); ax=axes(f);
imagesc(ax,G.graph); axis(ax,'image'); caxis(ax,[0 1]);
xlabel(ax,'Cause j'); ylabel(ax,'Effect i');
set(ax,'XTick',[1 5:5:r],'YTick',[1 5:5:r]);
colormap(f,[1 1 1; 0 0 0]); % 0: white (absent), 1: black (present). apply_taira_style(f);
set(ax,'FontName','Times New Roman','FontSize',26,'FontWeight','bold','LineWidth',1.3);
set([ax.XLabel ax.YLabel],'FontName','Times New Roman','FontSize',28,'FontWeight','bold');
end
