% Why does the tdm26 click fail at n=701? Throw or divergence?
L=load('refit_c1broad.mat'); pa=L.R.pa;
for nn=[1401 701]
    p=pa; p.n=nn;
    p.isv=unique(max(1,min(nn,round(((1391:-10:11)/1401)*nn))),'stable');
    fprintf('\nn=%4d  numel(isv)=%d  isv range [%d..%d]\n', nn, numel(p.isv), min(p.isv), max(p.isv));
    err='';
    t0=tic;
    try, evalc('S=tdm26(0,p,0,0);'); catch e, err=e.message; end
    el=toc(t0);
    if (~isempty(err))
        fprintf('  THREW after %.1f s: %s\n', el, err);
    else
        nf=sum(~isfinite(S.d1(:)));
        fprintf('  ran %.1f s | d1 %d/%d non-finite | max|d1| %.3e | ped max %.3e\n', ...
                el, nf, numel(S.d1), max(abs(S.d1(isfinite(S.d1)))), max(abs(S.ped)));
        fprintf('  size(d1)=%s  numel(S.f)=%d\n', mat2str(size(S.d1)), numel(S.f));
    end
end
disp('N701_DIAG_DONE');
