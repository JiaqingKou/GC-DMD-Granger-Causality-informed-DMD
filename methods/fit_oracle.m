function A=fit_oracle(data,lambda,A_true,opts)
r=size(A_true,1); W=double(A_true==0); W(1:r+1:end)=0;
A=optimize_operator('solve',data,lambda,W,opts);
end
