function D=compute_cavity_figure_data(D,N,indices,case_name)
%COMPUTE_CAVITY_FIGURE_DATA Convert precomputed streamfunction fields to vorticity.
% Recursive prediction, physical reconstruction and errors are computed in main.
nt=D.nt;
% Use exactly the original per-case interpolation grid and operator.
switch lower(case_name)
    case 'cavity13k', grid_n=301; use_sorted=false;
    case 'cavity19k', grid_n=401; use_sorted=true;
    case 'cavity30k', grid_n=201; use_sorted=true;
    otherwise, error('Unknown cavity figure case: %s',case_name);
end
Grid=CavityGridOperators(N); D2=Grid.D^2;
[xi,yi]=meshgrid(linspace(-1,1,grid_n));
D.xi=xi; D.yi=yi;
indices=unique(round(indices(:)'),'stable');
assert(all(indices>=1 & indices<=nt),'Figure snapshot index exceeds available data.');
D.indices=indices; D.fields=cell(1,nt);
for it=indices
    psi={D.Xtruth(:,it),D.Xd(:,it),D.Xg(:,it)};
    w=cell(1,3);
    for j=1:3
        P=reshape(real(psi{j}),N+1,N+1);
        W=-(D2*P+P*D2');
        if use_sorted
            [x,ord]=sort(Grid.x);
            w{j}=interp2(x,x,W(ord,ord),xi,yi,'spline');
        else
            [xx,yy]=meshgrid(Grid.x,Grid.x);
            w{j}=interp2(xx,yy,W,xi,yi,'spline');
        end
    end
    D.fields{it}={w{1},w{2},w{3},w{2}-w{1},w{3}-w{1}};
end
end
