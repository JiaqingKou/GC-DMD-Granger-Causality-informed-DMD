function M = modal_analysis_kou(lam,Phi,b,N,dt)
%MODAL_ANALYSIS_KOU Rank already-computed DMD modes by Kou contribution.
% Eigenvalues, eigenvectors, modal amplitudes and physical modes are formed
% explicitly in main; this function only ranks and packages the results.
phi2=sum(abs(Phi).^2,1).';
k=(0:N-1).';
I=zeros(numel(lam),1);
for j=1:numel(lam)
    I(j)=dt*sum(abs(b(j)*(lam(j).^k)))*phi2(j);
end
[Is,ord]=sort(real(I),'descend');
frac=Is/max(sum(Is),eps);
ls=lam(ord); bs=b(ord);
M.lambda=ls; M.integrated=Is; M.fraction=frac; M.order=ord;
M.frequency=angle(ls)/(2*pi*dt);
M.dt=dt;
M.b=bs; % complex initial modal amplitudes, sorted consistently
M.Phi=Phi(:,ord); % full-field DMD modes, sorted by contribution
M.effective_modes=1/max(sum(frac.^2),eps);
M.table=table((1:numel(lam)).',ord,real(ls),imag(ls),abs(ls),log(max(abs(ls),eps)),angle(ls)/(2*pi),M.frequency,abs(bs),Is,frac, ...
    'VariableNames',{'ModeRank','EigIndex','EigenvalueReal','EigenvalueImag','Magnitude','GrowthPerStep','FrequencyCyclesPerStep','Frequency','InitialAmplitude','IntegratedContribution','ContributionFraction'});
end
