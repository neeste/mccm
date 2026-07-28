% Did the fit DECOMPRESS the CF map, or just lower maperr some other way?
% Compare BF-per-place: native m=4 (0.32-2.10 kHz), fitted m=4, and m=2 (target).
ISV=[1136 1005 840 655 466 273 80];
L=load('parfit26_c4map.mat'); paf=L.R.pa;
cfg={ {'m=2 (target)', modpar26(2)}, {'m=4 native', modpar26(4)}, {'m=4 FITTED', paf} };
fprintf('\n  x/L:      0.811   0.717   0.600   0.468   0.333   0.195   0.057   | range(oct)\n');
for c=1:numel(cfg)
    nm=cfg{c}{1}; pa=cfg{c}{2}; pa.isv=ISV;
    evalc('S=tdm26(0,pa,0,0);');
    nf=numel(S.f); f=S.f(:); P=fft(S.ped); P=P(1:nf); bf=nan(1,numel(ISV));
    for i=1:numel(ISV)
        D=fft(S.d1(:,i)); H=D(1:nf)./max(abs(P),eps); a=abs(H).*f; a(~(f>0.15&f<18))=0;
        [pk,ip]=max(a); if(pk>0), bf(i)=f(ip); end
    end
    v=bf(isfinite(bf)&bf>0); rng=NaN; if(numel(v)>=2), rng=log2(max(v)/min(v)); end
    mono=''; if(numel(v)>=3 && ~all(diff(bf(isfinite(bf)))>0)), mono=' NONMONO'; end
    fprintf('  %-12s %s | %.2f%s\n', nm, num2str(bf,'%7.2f'), rng, mono);
end
% also report the fitted chsz and key slopes
fprintf('\n  fitted chsz = [%s]   (native [0.95 0.05 1.0 0.05])\n', num2str(paf.chsz,'%.3f '));
fprintf('  fitted k1e=%.3f m1e=%.3f  (native k1e=%.3f m1e=%.3f)\n', ...
        paf.k1e, paf.m1e, modpar26(4).k1e, modpar26(4).m1e);
disp('C4MAP_BF_DONE');
