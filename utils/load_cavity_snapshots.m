function D=load_cavity_snapshots(data_file,struct_name,field_name,n_model,pred_horizon)
%LOAD_CAVITY_SNAPSHOTS Load raw snapshots and grid/time metadata only.
S=load(data_file,struct_name);
assert(isfield(S,struct_name),'Missing struct %s.',struct_name);
VF=S.(struct_name);
assert(isfield(VF,field_name),'Missing field %s.%s.',struct_name,field_name);
Xraw=double(VF.(field_name)); nt=n_model+pred_horizon;
assert(size(Xraw,2)>=nt,'Need at least %d snapshots.',nt);
Xraw=Xraw(:,1:nt);
if isfield(VF,'N'), N=double(VF.N); else, N=round(sqrt(size(Xraw,1)))-1; end
assert((N+1)^2==size(Xraw,1),'Field size is incompatible with cavity grid N.');
if isfield(VF,'Time') && numel(VF.Time)>=nt
    time=reshape(double(VF.Time(1:nt)),1,[]);
    assert(all(isfinite(time)) && all(diff(time)>0),'Time must be finite and increasing.');
    dt=mean(diff(time));
else
    dt=0.1; time=(0:nt-1)*dt;
    warning('VelocityField.Time missing: using nondimensional elapsed time (k-1)*0.1.');
end
D=struct('Xraw',Xraw,'N',N,'time',time,'dt',dt);
fprintf('Field: %d x %d, N=%d, dt=%.6g\n',size(Xraw,1),nt,N,dt);
end
