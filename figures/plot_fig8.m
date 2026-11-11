function fig=plot_fig8(G,lambda_grid,val_gc,val_dmd,lambda_pred,Ad,Ag,cmax,gmax,outdir)
% One canvas: GC | lambda selection; DMD operator | GC-DMD operator.
C=taira_colors(); fig=figure('Name','GC selection and operators','Color','w','Position',[80 60 1100 1100]);
pp={[.10 .60 .30 .30],[.59 .60 .30 .30],[.10 .12 .30 .30],[.59 .12 .30 .30]};
axs=gobjects(1,4);
for k=1:4
    ax=axes(fig,'Position',pp{k}); axs(k)=ax;
    if k==2
        pos=lambda_grid>0;
        semilogx(ax,lambda_grid(pos),val_gc(pos),'-o','Color',C.red,'LineWidth',1.5,'MarkerSize',3); hold(ax,'on');
        yline(ax,val_dmd,'--','Color',C.gray,'LineWidth',1.1);
        if lambda_pred>0
            xline(ax,lambda_pred,'--','Color',C.blue,'HandleVisibility','off');
            yl=ylim(ax);
            text(ax,lambda_pred,val_dmd-.025*diff(yl),'  \lambda_{pred}', ...
                'Interpreter','tex','Rotation',90,'VerticalAlignment','top', ...
                'Color',C.blue,'FontSize',18,'Clipping','on');
        end
        lp=lambda_grid(pos);
        lo=floor(log10(min(lp))); hi=ceil(log10(max(lp)));
        xlim(ax,[1e-16 1e-3]);
        exponents=[-16 -12 -9 -6 -3];
        xticks(ax,10.^exponents); % approximately three decade labels
        ax.XTickLabelRotation=0;
        xlabel(ax,'\lambda'); ylabel(ax,'RMSE_{val}(\it A\rm)'); ylim(ax,[0 8e-4]); yticks(ax,0:2e-4:8e-4);
        legend(ax,'GC-DMD','DMD (\lambda=0)','Location','best');
    else
        if k==1
            imagesc(ax,G); caxis(ax,[0 gc_display_limit(G,gmax)]); colormap(ax,gc_white_red(256));
            xlabel(ax,'Cause \itj'); ylabel(ax,'Effect \iti'); tag='GC';
        else
            if k==3, A=Ad; else,A=Ag;end
            imagesc(ax,A); caxis(ax,[-cmax cmax]); colormap(ax,taira_diverging(256));
            xlabel(ax,'Source \itj'); ylabel(ax,'Target \iti'); tag='A_{ij}';
        end
        axis(ax,'image'); paper_matrix_ticks(ax,size(G,1));
        cb=colorbar(ax); cb.Position=[pp{k}(1)+.315 pp{k}(2) .012 .30]; cb.Ticks=linspace(ax.CLim(1),ax.CLim(2),5);

    end
end
apply_taira_style(fig,'cavity19k');
for k=1:4
    axs(k).PositionConstraint='innerposition'; axs(k).Position=pp{k}; pbaspect(axs(k),[1 1 1]);
    % A second xlabel line reserves space for the panel letter.
    axis_label=axs(k).XLabel.String;
    axs(k).XLabel.String={axis_label,sprintf('(%c)','a'+k-1)};
end
legs=findall(fig,'Type','legend'); for k=1:numel(legs),legs(k).FontSize=17;end
cbs=findall(fig,'Type','colorbar'); for k=1:numel(cbs),cbs(k).FontSize=17;end
% Rendering is deferred; no forced graphics synchronization.
for k=[1 3 4]
    axs(k).Position=pp{k};
end
set(fig,'PaperUnits','inches','PaperPositionMode','manual','PaperPosition',[0 0 8 7.2],'PaperSize',[8 7.2]);
if ~exist(outdir,'dir'),mkdir(outdir);end
% Automatic EPS export disabled; display first.
% No automatic PDF export.
end

function a=gc_display_limit(G,requested)
% Display only: omit self terms; use 99th percentile of positive entries.
if nargin>1 && ~isempty(requested), a=requested; else
    v=G(~eye(size(G))); v=v(isfinite(v)&v>0);
    if isempty(v), a=1; else, a=prctile(v,99); end
end
if ~isscalar(a)||~isfinite(a)||a<=0, error('GC color maximum must be positive.'); end
end

function cmap = gc_white_red(n)
%GC_WHITE_RED White at zero; increasingly saturated red for positive GC.
if nargin < 1 || isempty(n), n = 256; end
n = max(2,round(n));
t = linspace(0,1,n)';
% Match the operator-style white background and red active couplings.
anchor = [1 1 1; 0.98 0.89 0.88; 0.91 0.53 0.51; 0.70 0.11 0.12];
cmap = interp1(linspace(0,1,size(anchor,1)),anchor,t,'linear');
end

function paper_matrix_ticks(ax,n)
% At most six integer ticks, including both endpoints.
t=unique(round(linspace(1,n,min(n,6))));
set(ax,'XTick',t,'YTick',t,'XTickLabelRotation',0,'YTickLabelRotation',0);
end
