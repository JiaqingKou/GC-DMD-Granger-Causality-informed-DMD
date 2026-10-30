function finalize_paper_figures(figs)
% Display only: math notation and explicit endpoint ticks, including both yyaxis sides.
for f=reshape(figs,1,[])
 axs=findall(f,'Type','axes');
 for ax=reshape(axs,1,[])
  set(ax,'XGrid','off','YGrid','off','XMinorGrid','off','YMinorGrid','off', ...
   'XTickLabelRotation',0,'FontName','Times New Roman','FontWeight','bold');
  imgs=findall(ax,'Type','image');
  if ~isempty(imgs)
   sz=size(imgs(1).CData);nx=sz(2);ny=sz(1);
   set(ax,'XTick',unique([1 5:5:nx nx]),'YTick',unique([1 5:5:ny ny]));
   % Hidden repeated y labels remain hidden; visible labels receive the 1 tick.
  else
   % Ordinal plots have user-supplied tick labels; preserve those labels.
   if strcmp(ax.XTickLabelMode,'auto') || strcmp(ax.XScale,'log')
    v=endpoint_ticks(ax.XLim,ax.XTick,ax.XScale);
    ax.XLim=[v(1) v(end)];ax.XTick=v;
    if strcmp(ax.XScale,'log'),ax.XTickLabelMode='auto';end
   end
   for j=1:numel(ax.YAxis)
    yr=ax.YAxis(j);v=endpoint_ticks(yr.Limits,yr.TickValues,yr.Scale);
    yr.Limits=[v(1) v(end)];yr.TickValues=v;
   end
  end
  math_text(ax.XLabel);
  for j=1:numel(ax.YAxis),math_text(ax.YAxis(j).Label);end
  math_text(ax.Title);
 end
 cbs=findall(f,'Type','colorbar');
 for cb=reshape(cbs,1,[]),math_text(cb.Label);end
 txt=findall(f,'Type','text');
 for h=reshape(txt,1,[]),math_text(h);end
 legs=findall(f,'Type','legend');
 for lg=reshape(legs,1,[])
  ss=lg.String;
  for j=1:numel(ss),ss{j}=math_string(ss{j});end
  lg.String=ss;
  if any(cellfun(@(s)contains(s,'\bf'),ss))
   for j=1:numel(ss)
    if strcmp(ss{j},'$F1$'),ss{j}='\it F\rm\bf 1';end
   end
   lg.String=ss;lg.Interpreter='tex';lg.FontWeight='bold';
  elseif any(cellfun(@(s)contains(s,'$'),ss)),lg.Interpreter='latex';
  else,lg.Interpreter='tex';lg.FontWeight='bold';end
 end
end
end
function v=endpoint_ticks(lim,v,scale)
% Uniform ticks: expand to round limits rather than append a short last interval.
if strcmp(scale,'log')
 lo=floor(log10(lim(1))+1e-10);hi=ceil(log10(lim(2))-1e-10);
 stride=max(1,ceil((hi-lo)/7));
 ex=lo:stride:hi;
 if ex(end)<hi,ex=linspace(lo,hi,ceil((hi-lo)/stride)+1);end
 v=10.^ex;
else
 d=diff(lim)/5;base=10^floor(log10(max(d,realmin)));
 choices=[1 2 2.5 5 10]*base;
 step=choices(find(choices>=d*(1-1e-10),1));
 lo=floor(lim(1)/step+1e-9)*step;hi=ceil(lim(2)/step-1e-9)*step;
 v=lo:step:hi;v(abs(v)<step*1e-9)=0;
end
end
function math_text(h)
s=h.String;
if iscell(s)
 for j=1:numel(s),s{j}=math_string(s{j});end
else
 s=math_string(s);
end
h.String=s;
if ischar(s)&&contains(s,'\bf'),h.Interpreter='tex';
 if contains(s,'\lambda'),h.FontWeight='normal';end
else,h.Interpreter='latex';end
end
function s=math_string(s)
if ~ischar(s)||isempty(s)||contains(s,'$'),return;end
switch s
 case 'Source j',s='\bfSource \itj\rm\bf';
 case 'Target i',s='\bfTarget \iti\rm\bf';
 case 'Cause j',s='\bfCause \itj\rm\bf';
 case 'Effect i',s='\bfEffect \iti\rm\bf';
 case '\lambda',s='$\lambda$';
 case 'i',s='$i$';
 case 'j',s='$j$';
 case 'p',s='$p$';
 case 'N',s='$N$';
 case 'k',s='$k$';
 case '|A_{ij}^{true}|',s='$|A_{ij}^{\mathrm{true}}|$';
 case 'Gtilde_{ij}',s='$\widetilde{G}_{ij}$';
 case 'RMSEval(A)',s='\bfRMSE_{val}(\itA\rm\bf)';
 case 'F1',s='\itF\rm\bf1';
 case 'EA',s='$\|\mathbf{A}-\mathbf{A}^{\mathrm{true}}\|_F/\|\mathbf{A}^{\mathrm{true}}\|_F$';
 case 'Epred',s='$\|\widehat{\mathbf{Q}}-\mathbf{Q}\|_F/\|\mathbf{Q}\|_F$';
 case 'E(k)',s='$E(k)$';
 case 'density',s='$\mathrm{nnz}_{i\ne j}(\mathbf{A})/[r(r-1)]$';
 case 'L1off',s='$\sum_{i\ne j}|A_{ij}|$';
 case 'rho(A)',s='$\rho(\mathbf{A})$';
 case 'A_{true}',s='$\mathbf{A}^{\mathrm{true}}$';
 case 'A_{GC-DMD,ij} - A_{true,ij}',s='$A_{ij}^{\mathrm{GC-DMD}}-A_{ij}^{\mathrm{true}}$';
 otherwise
  if contains(s,'\lambda')
   % Native TeX preserves Times New Roman bold for annotations and numbers.
   for tag={'pred','gc','weak','strong','over'}
    s=strrep(s,['\lambda_{' tag{1} '}'],['\lambda_{\bf ' tag{1} '}']);
   end
   s=strrep(s,'=','\rm\bf =');
  end
end
end
