function [A,info]=fit_l1_dmd(data,lambda,opts)
%FIT_L1_DMD Uniform-L1 DMD, using the same FISTA optimizer as GC-DMD.
% Appendix A penalizes all coefficients, including diagonal entries.
r=size(data,2);
W=ones(r);
[A,info]=optimize_operator('solve',data,lambda,W,opts);
end
