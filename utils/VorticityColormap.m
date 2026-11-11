function cmap=VorticityColormap(N)
% Original Arbabi/Mezic color values, without creating or closing any figure.
if nargin<1, N=64; end
assert(mod(N,16)==0,'Colormap size must be divisible by 16.');
cmap=jet(N);
for i=N/8+1:7*N/16
    cmap(i,:)=[(i-N/8)/(5*N/16),(i-N/8)/(5*N/16),1];
end
cmap(7*N/16+1:9*N/16,:)=1;
for i=9*N/16+1:7*N/8
    cmap(i,:)=[1,1-(i-9*N/16)/(5*N/16),1-(i-9*N/16)/(5*N/16)];
end
end
