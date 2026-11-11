function [Grid,varargout] = CavityGridOperators(N)
%CAVITYGRIDOPERATORS Self-contained Chebyshev cavity grid and operators.
% No external paths or third-party helper functions are required.

Grid = CollocationGrid_q(N);

if nargout == 2
    varargout{1} = CreateOperators_psi(Grid.D,eye(N+1));
end

[~,Grid.wc] = clencurt_local(N);
Grid.W = kron(Grid.wc,Grid.wc);
end

function Grid = CollocationGrid_q(N)
[Grid.D,Grid.x] = cheb_local(N);
[xx,yy] = meshgrid(Grid.x,Grid.x);
Grid.xx = xx(:);
Grid.yy = yy(:);

Grid.Walls      = find(abs(xx)==1 | abs(yy)==1);
Grid.LeftWall   = find(xx==-1);
Grid.RightWall  = find(xx==1);
Grid.BottomWall = find(yy==-1);
Grid.TopWall    = find(yy==1);

I = eye(N+1);
Grid.D2 = diag(1-Grid.x.^2)*(Grid.D^2) - 4*diag(Grid.x)*Grid.D - 2*I;
Grid.D4 = diag(1-Grid.x.^2)*(Grid.D^4) - 8*diag(Grid.x)*(Grid.D^3) - 12*(Grid.D^2);
end

function [D,x] = cheb_local(N)
%CHEB_LOCAL Chebyshev differentiation matrix on N+1 Lobatto points.
if N == 0
    D = 0;
    x = 1;
    return
end
x = cos(pi*(0:N)'/N);
c = [2; ones(N-1,1); 2] .* ((-1).^(0:N)');
X = repmat(x,1,N+1);
dX = X - X.';
D = (c*(1./c).') ./ (dX + eye(N+1));
D = D - diag(sum(D,2));
end

function [x,w] = clencurt_local(N)
%CLENCURT_LOCAL Clenshaw-Curtis quadrature nodes and weights.
if N == 0
    x = 1;
    w = 2;
    return
end
theta = pi*(0:N)'/N;
x = cos(theta);
w = zeros(1,N+1);
ii = 2:N;
v = ones(N-1,1);
if mod(N,2) == 0
    w(1) = 1/(N^2-1);
    w(N+1) = w(1);
    for k = 1:N/2-1
        v = v - 2*cos(2*k*theta(ii))/(4*k^2-1);
    end
    v = v - cos(N*theta(ii))/(N^2-1);
else
    w(1) = 1/N^2;
    w(N+1) = w(1);
    for k = 1:(N-1)/2
        v = v - 2*cos(2*k*theta(ii))/(4*k^2-1);
    end
end
w(ii) = 2*v/N;
end

function Operators_psi = CreateOperators_psi(D,I)
%CREATEOPERATORS_PSI Derivative operators for vectorized streamfunction.
Operators_psi.Dx = kron(D,I);
Operators_psi.Dy = kron(I,D);
D2 = D^2;
D4 = D2^2;
Operators_psi.del2 = kron(D2,I) + kron(I,D2);
Operators_psi.del4 = kron(D4,I) + kron(I,D4) + 2*kron(D2,I)*kron(I,D2);
end
