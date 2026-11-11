function [W,GC_raw,Sig,C]=build_gc_weights(input,settings)
%BUILD_GC_WEIGHTS GC-DMD penalty: W=1-C, with unpenalized diagonal.
% An optional settings.input_type='data' estimates GC for a refit.
if nargin<2, settings=struct; end
GC_raw=[]; Sig=[];
if isfield(settings,'input_type') && strcmp(settings.input_type,'data')
    [GC_raw,~,Sig,C]=estimate_gc_var1(input,'alpha',settings.alpha, ...
        'fdr_q',settings.alpha,'norm_mode',settings.norm_mode);
else
    C=input;
end
r=size(C,1);
W=1-C;
W(1:r+1:end)=0;
end
