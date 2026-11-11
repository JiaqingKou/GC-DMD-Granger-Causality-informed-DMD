function [U,s,V]=compute_svd(X)
%COMPUTE_SVD Economy snapshot SVD of the supplied matrix (never centers X).
G=X'*X;
[V,d]=eig((G+G')/2,'vector');
[s,idx]=sort(sqrt(max(real(d),0)),'descend');
V=V(:,idx);
keep=s>max(size(X))*eps(max(s));
s=s(keep); V=V(:,keep);
U=(X*V)./s';
end
