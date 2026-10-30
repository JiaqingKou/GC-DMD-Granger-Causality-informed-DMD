function [W,GC_raw,Sig,C] = build_gc_weights(input,settings)
%BUILD_GC_WEIGHTS Construct GC-informed operator penalty weights.
%   W = build_gc_weights(C,settings) uses an existing GC weight matrix C.
%   [W,GC_raw,Sig,C] = build_gc_weights(data,settings) estimates GC from
%   snapshots (N x r), then constructs the penalty. In both cases the
%   diagonal is unpenalized. Set settings.input_type to ''gc'' or ''data''
%   explicitly to avoid ambiguity for square input matrices.
%
if ~isfield(settings,'eta'), error('settings.eta is required.'); end
if isfield(settings,'input_type')
    input_type=lower(settings.input_type);
else
    if size(input,1)==size(input,2)
        input_type='gc';
    else
        input_type='data';
    end
end
GC_raw=[]; Sig=[];
switch input_type
    case 'gc'
        C=input;
    case 'data'
        [GC_raw,~,Sig,C,~]=estimate_gc_var1(input,'alpha',settings.alpha, ...
            'multiple_testing','bh-fdr','fdr_q',settings.alpha, ...
            'norm_mode',settings.norm_mode, ...
            'intercept',get_gc_intercept(settings));
    otherwise
        error('settings.input_type must be ''gc'' or ''data''.');
end
r=size(C,1);
if size(C,2)~=r, error('GC matrix must be square.'); end
W=1-settings.eta*C;
W(1:r+1:end)=0;
end

function tf = get_gc_intercept(settings)
% Keep GC model specification identical during training and final refit.
tf = true;
if isfield(settings,'intercept'), tf = settings.intercept; end
end
