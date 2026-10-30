function A=fit_gc_dmd(data,lambda,C,opts)
r=size(C,1); W=1-C; W(1:r+1:end)=0;
A=optimize_operator('solve',data,lambda,W,opts);
end
