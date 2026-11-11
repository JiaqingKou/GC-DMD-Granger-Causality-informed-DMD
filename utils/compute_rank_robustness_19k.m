function D=compute_rank_robustness_19k(Xraw,mu,Xc,V,sig,cum_energy,n_model,n_fit,n_val,pred_horizon,n_lambda,fista_opts,energy_targets)
% Original Figure 12 numerical rank sweep, moved unchanged from plot_fig12.
if nargin<13 || isempty(energy_targets)
    energy_targets=[0.99 0.999 0.9999 0.99999];
end
C=taira_colors();
nt=min(size(Xraw,2),n_model+pred_horizon);
H=nt-n_model;
assert(H>=1,'No held-out snapshots are available for rank robustness.');

% Convert energy targets to unique ranks, preserving increasing order.
ranks=zeros(size(energy_targets));
for k=1:numel(energy_targets)
    j=find(cum_energy>=energy_targets(k),1,'first');
    if isempty(j), j=numel(sig); end
    ranks(k)=j;
end
[ranks,ia]=unique(ranks,'stable');
energy_targets=energy_targets(ia);
nr=numel(ranks);

rhoD=zeros(nr,1); rhoG=zeros(nr,1);
meanD=zeros(nr,1); meanG=zeros(nr,1);
endD=zeros(nr,1); endG=zeros(nr,1);
lambdaGC=zeros(nr,1);
retained=zeros(nr,1);
densityD=zeros(nr,1); densityG=zeros(nr,1);
sparsity_tol_rel=1e-3;

fprintf('\n=========== POD-rank / long-prediction robustness ===========\n');
for kk=1:nr
    rr=ranks(kk); retained(kk)=cum_energy(rr);
    fprintf('  [%d/%d] target=%.5f, r=%d, retained=%.8f\n', ...
        kk,nr,energy_targets(kk),rr,retained(kk));

    Urr=Xc(:,1:n_model)*(V(:,1:rr)./sig(1:rr)');
    [Urr,~]=qr(Urr,0);
    Qr=(Urr'*Xc(:,1:nt))';
    Qm=Qr(1:n_model,:);
    tr=Qm(1:n_fit,:);
    va=Qm(n_fit+1:n_fit+n_val,:);

    lg=build_lambda_grid(Qm(1:n_fit,:),n_lambda,'cavity19k');
    R=identify_operators(Qm,n_fit,lg,fista_opts);
    Ad=R.A_dmd; Ag=R.A_gc; lambdaGC(kk)=R.lambda_gc;

    rhoD(kk)=max(abs(eig(Ad)));
    rhoG(kk)=max(abs(eig(Ag)));
    sd=operator_sparsity_stats(Ad,sparsity_tol_rel);
    sg=operator_sparsity_stats(Ag,sparsity_tol_rel);
    densityD(kk)=sd.density; densityG(kk)=sg.density;

    % Recursive HELD-OUT prediction. Start from the measured state at n_model.
    qd=Qr(n_model,:)'; qg=qd;
    errD=zeros(H,1); errG=zeros(H,1);
    for h=1:H
        qd=Ad*qd; qg=Ag*qg;
        it=n_model+h;
        xd=mu+Urr*qd; xg=mu+Urr*qg;
        xt=Xraw(:,it);
        den=norm(xt-mu); if den<eps, den=max(norm(xt),eps); end
        errD(h)=norm(xd-xt)/den;
        errG(h)=norm(xg-xt)/den;
    end
    meanD(kk)=100*mean(errD); meanG(kk)=100*mean(errG);
    endD(kk)=100*errD(end); endG(kk)=100*errG(end);
    fprintf('       rho: DMD %.6f, GC %.6f | mean pred err: DMD %.4e, GC %.4e\n', ...
        rhoD(kk),rhoG(kk),meanD(kk),meanG(kk));
end

D=struct('ranks',ranks,'rhoD',rhoD,'rhoG',rhoG,'meanD',meanD,'meanG',meanG,'H',H);
end
