function cmap = supp_taira_diverging(n)
%TAIRA_DIVERGING Blue-white-red diverging map with white at zero.
if nargin<1, n=257; end
n1=floor(n/2); n2=n-n1;
blue=[0.12 0.35 0.70]; white=[1 1 1]; red=[0.75 0.18 0.16];
a=linspace(0,1,n1+1)'; b=linspace(0,1,n2+1)';
left=blue.*(1-a)+white.*a;
right=white.*(1-b)+red.*b;
cmap=[left(1:end-1,:); right];
end
