function [fig9a,fig9b,fig9c]=plot_fig9(Md,Mg,dt,models,N,n_show)
% Figure 9 (a,b,c): source drawing statements preserved from original 19k project.
C=taira_colors(); if nargin<3 || isempty(dt),dt=1;end

fig=figure('Name','Ritz eigenvalues','Color','w','Position',[40 80 620 520]);
th=linspace(0,2*pi,500);
ax=axes(fig); plot(ax,cos(th),sin(th),'k--','LineWidth',1.1); hold(ax,'on');
plot(ax,real(Md.lambda),imag(Md.lambda),'o','Color',C.blue,'MarkerFaceColor','w','LineWidth',1.2,'MarkerSize',6);
plot(ax,real(Mg.lambda),imag(Mg.lambda),'x','Color',C.red,'LineWidth',1.5,'MarkerSize',7); axis(ax,'square'); xlim(ax,[0 1.2]); xticks(ax,0:.4:1.2); ylim(ax,[-.8 .8]); yticks(ax,-.8:.4:.8); xlabel(ax,'Re(\lambda)'); ylabel(ax,'Im(\lambda)'); title(ax,'(a) Ritz eigenvalues'); legend(ax,'Unit circle','DMD','GC-DMD','Location','best');
wd=log(Md.lambda)/dt; wg=log(Mg.lambda)/dt; fd=abs(imag(wd))/(2*pi); fg=abs(imag(wg))/(2*pi); gd=real(wd); gg=real(wg);
align_paper_panel(fig,'cavity19k'); apply_taira_style(fig,'cavity19k');

fig9a=fig;
fd2=angle(Md.lambda)/(2*pi*dt); fg2=angle(Mg.lambda)/(2*pi*dt); kd=fd2>=-1e-12; kg=fg2>=-1e-12;
fig=figure('Name','Frequency contribution spectrum','Color','w','Position',[120 80 620 520]); ax=axes(fig); stem(ax,fd2(kd),Md.fraction(kd),'o','Color',C.blue,'MarkerFaceColor','w','LineWidth',1.25,'MarkerSize',4); hold(ax,'on'); stem(ax,fg2(kg),Mg.fraction(kg),'s','Color',C.red,'MarkerFaceColor','w','LineWidth',1.25,'MarkerSize',4); xlim(ax,[0 1]); xticks(ax,0:.25:1); yl=ylim(ax); ymax=0.12; ylim(ax,[0 ymax]); yticks(ax,0:0.04:ymax); xlabel(ax,'\itf'); ylabel(ax,'\eta'); title(ax,'(c) Frequency spectrum'); legend(ax,'DMD','GC-DMD','Location','best');
align_paper_panel(fig,'cavity19k'); apply_taira_style(fig,'cavity19k');
fig9b=fig;
original_modal_atlas_19k(models,N,n_show,'Contribution-ranked vorticity modal atlas');
fig9c=gcf;
end

function stats = original_modal_atlas_19k(models,N,n_show,ttl)
%PLOT_CAVITY_MODAL_ATLAS Contribution-ranked one-conjugate-side vorticity modes.
% Combined atlas: each method occupies one row, one conjugate side only.
if nargin<3 || isempty(n_show), n_show=9; end
if nargin<4, ttl='Contribution-ranked vorticity modal atlas'; end

[Grid,Operators]=CavityGridOperators(N); cmap=VorticityColormap(64);
[xx2,yy2]=meshgrid(Grid.x,Grid.x); [xi,yi]=meshgrid(linspace(-1,1,401));
stats = repmat(struct('name','','rank',[],'frequency',[], ...
    'contribution',[],'roughness',[]), numel(models), 1);
fig=figure('Name',ttl,'Color','w','Position',[70 70 1380 680]);
tl=tiledlayout(fig,numel(models),n_show,'TileSpacing','compact','Padding','compact');
add_vertical_method_labels(fig,{models.name});
tl.Units='normalized'; tl.OuterPosition=[0.085 0.02 0.91 0.96];
for im=1:numel(models)
    M=models(im).M; f=M.table.FrequencyCyclesPerStep;
    keep=find(f>=-1e-12); keep=keep(1:min(n_show,numel(keep)));
    rough=zeros(numel(keep),1); orig_rank=rough; freq=rough; contrib=rough;
    for k=1:numel(keep)
        rk=keep(k); ax=nexttile(tl,(im-1)*n_show+k);
        SF=real(M.Phi(:,rk)); w=-Operators.del2*SF; ww=reshape(real(w),N+1,N+1);
        www=interp2(xx2,yy2,ww,xi,yi,'spline'); [gy,gx]=gradient(ww);
        rough(k)=sum(abs(gx(:)).^2+abs(gy(:)).^2)/max(sum(abs(ww(:)).^2),eps);
        orig_rank(k)=rk; freq(k)=f(rk); contrib(k)=M.fraction(rk);
        a=prctile(abs(www(:)),99); if ~isfinite(a)||a<=0,a=max(abs(www(:)));end; if ~isfinite(a)||a<=0,a=1;end
        levels=linspace(-a,a,101); www=max(min(www,a),-a);
        draw_reference_flow(ax,xi,yi,www,a); axis(ax,'square'); box(ax,'on');
        colormap(ax,cmap); caxis(ax,[-a a]); xlim(ax,[-1 1]); ylim(ax,[-1 1]);
        xticks(ax,[-1 0 1]);yticks(ax,[-1 0 1]);xlabel(ax,'\itx'); ylabel(ax,'\ity');
        text(ax,0.5,1.04,sprintf('Mode %d',rk),'Units','normalized','HorizontalAlignment','center','FontName','Times New Roman','FontSize',15,'FontWeight','bold','Clipping','off');
    end
    stats(im)=struct('name',models(im).name,'rank',orig_rank,'frequency',freq,'contribution',contrib,'roughness',rough);
end
apply_taira_style(fig,'cavity19k');
Tout=table();
for im=1:numel(models)
    ni=numel(stats(im).rank);
    Ti=table(repmat(string(stats(im).name),ni,1),stats(im).rank,stats(im).frequency,stats(im).contribution,stats(im).roughness, ...
        'VariableNames',{'Method','ContributionRank','Frequency','Contribution','SpatialRoughness'});
    Tout=[Tout;Ti]; %#ok<AGROW>
end
% No CSV output: paper figure only.
end

function add_vertical_method_labels(fig,names)
% Dedicated in-canvas strip: vertical labels do not extend beyond the figure.
ax=axes('Parent',fig,'Units','normalized','Position',[0.055 0.02 0.035 0.96], ...
    'Visible','off','XLim',[0 1],'YLim',[0 1],'Tag','PaperRowLabels','HandleVisibility','off');
for row=1:numel(names)
    text(ax,0.5,1-(row-0.5)/numel(names),names{row}, ...
        'Rotation',90,'HorizontalAlignment','center','VerticalAlignment','middle', ...
        'FontName','Times New Roman','FontSize',24,'FontWeight','bold','Clipping','on');
end
end

