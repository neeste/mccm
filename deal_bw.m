function b=deal_bw(A,f,ip,thr)
nf=numel(A); lo=ip; while (lo>1 && A(lo)>thr), lo=lo-1; end
hi=ip; while (hi<nf && A(hi)>thr), hi=hi+1; end
b=f(hi)-f(lo);
