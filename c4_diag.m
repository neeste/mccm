% Why does the fitted heavy-m2 model give NaN in the tdm26 click? Divergence, or
% a degenerate/folded map? fdm26 maperr said 464 (decent) but tdm26 click failed.
L=load('parfit26_c4map_m2.mat'); paf=L.R.pa;
paf.isv=[1136 1005 840 655 466 273 80];
fprintf('fitted chsz raw = [%s]  (normalized sum=2: [%s])\n', ...
        num2str(paf.chsz,'%.3f '), num2str(paf.chsz*2/sum(paf.chsz),'%.3f '));
err='';
try, evalc('S=tdm26(0,paf,0,0);'); catch e, err=e.message; end
if(~isempty(err)), fprintf('tdm26 THREW: %s\n',err); disp('C4_DIAG_DONE'); return; end
nfin=sum(~isfinite(S.d1(:)));
vv={'finite','DIVERGES'}; fprintf('d1: %d/%d non-finite  ->  %s\n', nfin, numel(S.d1), vv{(nfin>0)+1});
nf=numel(S.f); f=S.f(:); P=fft(S.ped); P=P(1:nf);
fprintf('\n place  x/L    BF(kHz)   peakmag   in-band?\n');
for i=1:7
    D=fft(S.d1(:,i)); H=D(1:nf)./max(abs(P),eps); a=abs(H).*f;
    [pkall,ipa]=max(a);                       % peak over ALL freq (unfiltered)
    ab=a; ab(~(f>0.15&f<18))=0; [pk,ip]=max(ab);
    bf=NaN; if(pk>0), bf=f(ip); end
    fprintf('  %d   %.3f   %7.3f   %.2e   bf@%.2f (all-band peak @%.2f kHz)\n', ...
            i, paf.isv(i)/paf.n, bf, pk, f(ip), f(ipa));
end
% compare fdm26's view of the same model (the thing the fit optimized)
Rc=fdm26(struct('cfmap',1,'pa',paf,'flst',500*2.^(-1:5)));
fprintf('\n fdm26 peak place x/L per freq (250..16000 Hz):\n  %s\n', num2str(Rc.xpk_bm,'%6.3f '));
fprintf(' fdm26 maperr = %.1f\n', L.R.Rf.maperr);
disp('C4_DIAG_DONE');
