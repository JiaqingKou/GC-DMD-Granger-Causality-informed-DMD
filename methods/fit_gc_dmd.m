function A=fit_gc_dmd(data,lambda,C,opts)
W=build_gc_weights(C);
A=optimize_operator('solve',data,lambda,W,opts);
end
