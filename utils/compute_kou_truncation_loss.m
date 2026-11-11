function T=compute_kou_truncation_loss(Md,Mg,Ur,Xeval,max_modes,sampling_indices)
% Original Figure 11 numerical loss computation, without plotting.
if nargin<5 || isempty(max_modes),max_modes=115;end
n=size(Xeval,2);
if nargin<6 || isempty(sampling_indices),sampling_indices=1:n;end
assert(numel(sampling_indices)==n && all(sampling_indices>=1),'Invalid snapshot indices.');
Q=Ur'*Xeval;
assert(n>=1 && size(Ur,1)==size(Xeval,1),'Incompatible evaluation data.');
assert(norm(Xeval,'fro')>0,'Evaluation field has zero norm.');
mmax=min([max_modes,numel(Md.lambda),numel(Mg.lambda)]);
assert(mmax>=1,'No available DMD modes.');
% Orthogonal decomposition evaluates the full-field norm without repeatedly
% constructing nSpace-by-nSnapshots reconstructions for every truncation.
Xperp=Xeval-Ur*Q; floor2=norm(Xperp,'fro')^2; clear Xperp;
den=norm(Xeval,'fro'); powers=reshape(sampling_indices,1,[])-1; loss=zeros(mmax,2);
mods={Md,Mg};
for method=1:2
    M=mods{method}; Phir=Ur'*M.Phi(:,1:mmax);
    Qrec=complex(zeros(size(Q)));
    for m=1:mmax
        Qrec=Qrec+Phir(:,m)*(M.b(m)*(M.lambda(m).^powers));
        % Complex Frobenius norm follows Eq.(36), including an incomplete
        % conjugate pair at odd truncation counts. Do not discard its imaginary part.
        loss(m,method)=100*sqrt(norm(Q-Qrec,'fro')^2+floor2)/den;
    end
end
T=table((1:mmax)',loss(:,1),loss(:,2), ...
    'VariableNames',{'RetainedDMDModes','DMDLossPercent','GCDMDLossPercent'});
end
