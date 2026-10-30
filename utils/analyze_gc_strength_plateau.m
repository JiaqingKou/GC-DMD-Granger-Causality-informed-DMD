function R=analyze_gc_strength_plateau(GC_raw,Sig,varargin)
% Tail plateau: first point entering a near-zero band and remaining there.
% This is an explicitly declared heuristic, not a significance test.
p=inputParser; addParameter(p,'tail_level',0.01); addParameter(p,'tail_window',15); parse(p,varargin{:});
r=size(GC_raw,1); off=~eye(r); vals=GC_raw(off); scale=max(vals);
if scale<=0, error('All off-diagonal GC strengths are zero.'); end
[sorted,order]=sort(vals/scale,'descend'); n=numel(sorted);
% All later points must remain below the chosen near-zero band (sorted series).
knee=find(sorted<=p.Results.tail_level,1,'first');
if isempty(knee), knee=n; end
% If the tail has no sustained length, report rather than invent a plateau.
min_tail=min(p.Results.tail_window,n);
if n-knee+1<min_tail
 warning('GC:NoSustainedPlateau','No sustained near-zero tail at selected tolerance.');
end
% Place threshold at boundary midpoint to retain the last value above band.
if knee>1, tau_norm=(sorted(knee-1)+sorted(knee))/2;
else, tau_norm=sorted(1); end
tau=tau_norm*scale; graph=(GC_raw>tau)&off;
R.sorted_normalized=sorted; R.order=order; R.tau_raw=tau;
R.tau_normalized=tau_norm; R.knee=knee; R.method='near-zero tail';
R.graph=graph; R.scale=scale; R.sig_bh=logical(Sig)&off;
R.agreement=nnz(graph(off)==R.sig_bh(off))/nnz(off);
R.tail_level=p.Results.tail_level; R.tail_window=min_tail;
end
