function [GC,p,Sig,C,W,info]=infer_gc_weights(Q)
%INFER_GC_WEIGHTS Fixed paper GC settings: VAR(1), intercept, BH-FDR, rowsum.
[GC,p,Sig,C,info]=estimate_gc_var1(Q, ...
    'alpha',0.05,'multiple_testing','bh-fdr','fdr_q',0.05, ...
    'norm_mode','rowsum');
W=build_gc_weights(C);
end
