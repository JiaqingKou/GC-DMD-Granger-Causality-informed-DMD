function fig=plot_fig6(Ad,Ag,Md,Mg,dt,cmax,outdir) %#ok<INUSD>
% Ordinary axes allow explicit square-box alignment without tiledlayout conflicts.
C=taira_colors();
fig=figure('Name','Combined operator and spectral paper figure','Color','w', ...
    'Position',[60 60 1100 1100]);
% Positions are assigned by the common layout function after plotting.
axs=gobjects(1,4); ops={Ad,Ag};
for k=1:2
    ax=axes('Parent',fig,'Units','normalized');axs(k)=ax;
    imagesc(ax,ops{k});axis(ax,'image');caxis(ax,[-cmax cmax]);
    colormap(ax,taira_diverging(256));
    xlabel(ax,{'Source \itj\rm',sprintf('(%c)',char(96+k))},'Interpreter','tex');ylabel(ax,'Target \iti\rm','Interpreter','tex');
    xticks(ax,1:size(ops{k},2));yticks(ax,1:size(ops{k},1));
    cb=colorbar(ax,'eastoutside');setappdata(ax,'OperatorColorbar',cb);cb.Ticks=linspace(-cmax,cmax,5);
    cb.Label.String='';
    % The layout function keeps this colorbar adjacent to its image.
end
ax=axes('Parent',fig,'Units','normalized');axs(3)=ax;
th=linspace(0,2*pi,500);hc=plot(ax,cos(th),sin(th),'k--','LineWidth',1.1);hold(ax,'on');
hd=plot(ax,real(Md.lambda),imag(Md.lambda),'o','Color',C.blue,'MarkerFaceColor','w','LineWidth',1.2,'MarkerSize',6);
hg=plot(ax,real(Mg.lambda),imag(Mg.lambda),'x','Color',C.red,'LineWidth',1.5,'MarkerSize',7);
xlim(ax,[0.6 1.2]);ylim(ax,[-0.4 0.4]);
xticks(ax,0.6:0.2:1.2);yticks(ax,-0.4:0.2:0.4);
xlabel(ax,{'Re(\lambda)','(c)'});ylabel(ax,'Im(\lambda)');
legend(ax,[hc hd hg],{'Unit circle','DMD','GC-DMD'},'Location','southwest');
ax=axes('Parent',fig,'Units','normalized');axs(4)=ax;
fd=angle(Md.lambda)/(2*pi*dt);fg=angle(Mg.lambda)/(2*pi*dt);
kd=fd>=-1e-12;kg=fg>=-1e-12;
hd=stem(ax,fd(kd),Md.fraction(kd),'o','Color',C.blue,'MarkerFaceColor','w','LineWidth',1.4,'MarkerSize',5);hold(ax,'on');
hg=stem(ax,fg(kg),Mg.fraction(kg),'s','Color',C.red,'MarkerFaceColor','w','LineWidth',1.4,'MarkerSize',5);
fmax=max([0.5;fd(kd);fg(kg)]);fmax=ceil(fmax*10)/10;
xlim(ax,[0 fmax]);xticks(ax,linspace(0,fmax,3));
ymax=max([Md.fraction(kd);Mg.fraction(kg);eps]);
ytop=max(0.1,ceil(ymax*5)/5); ylim(ax,[0 ytop]); yticks(ax,linspace(0,ytop,4));
xlabel(ax,{'\itf\rm','(d)'},'Interpreter','tex');ylabel(ax,'\it\eta','Interpreter','tex');
legend(ax,[hd hg],{'DMD','GC-DMD'},'Location','northeast');
apply_taira_style(fig,'cavity13k');
for k=1:4
    axs(k).PositionConstraint='innerposition';
    if isprop(axs(k),'Toolbar'),axs(k).Toolbar.Visible='off';end
end
align_operator_spectral_boxes(fig,axs);
% Fixed layout: no resize callback or forced graphics synchronization.
set(fig,'SizeChangedFcn',[]);
% Display only: no automatic EPS/PDF export.
end

function align_operator_spectral_boxes(fig,axs)
% Recompute from fixed quadrant centers, never from previously resized boxes.
if ~isgraphics(fig) || any(~isgraphics(axs)),return;end
% Original v8 normalized plot boxes and spacing, on a square canvas.
w=0.32;h=0.32;
xc=[0.25 0.71];yc=[0.75 0.27];
for col=1:2
    for row=1:2
        ax=axs(col+2*(row-1));ax.Units='normalized';
        ax.PositionConstraint='innerposition';
        ax.Position=[xc(col)-w/2,yc(row)-h/2,w,h];
        pbaspect(ax,[1 1 1]);
    end
    at=axs(col);square=at.Position;
    cb=getappdata(at,'OperatorColorbar');
    if ~isempty(cb) && isgraphics(cb)
        cb.Units='normalized';cb.Position=[square(1)+w+0.002 square(2) 0.012 h];
        at.Position=square;
        lab=getappdata(at,'OperatorColorbarHeading');
        if ~isempty(lab) && isgraphics(lab)
            lab.Position=[cb.Position(1)-0.034 square(2)+h+0.006 0.08 0.03];
        end
    end
end
end
