% Is the +116 dB "amplifier" real, or a collapsed PASSIVE denominator? And is
% the fold=6.36 oct map genuinely broken? amp = peak@2k(active) - peak@2k(passive),
% so a tiny passive response inflates the difference without the active response
% being strong.
L=load('parfit26_c3amp.mat'); paf=L.R.pa;
ISV=1391:-10:11; XLO=0.05; XHI=0.85;
for tag={'native',modpar26(3); 'FITTED',paf}'
end
cfg={ {'native m=3', modpar26(3)}, {'FITTED c3amp', paf} };
for c=1:numel(cfg)
    nm=cfg{c}{1}; pa=cfg{c}{2}; pa.isv=ISV;
    fprintf('\n===== %s =====\n', nm);
    for og=[1 0]
        p=pa; if (p.m>=4), p.ohcgain=og; else, p.gam=og*pa.gam; end
        evalc('S=tdm26(0,p,0,0);');
        nf=numel(S.f); f=S.f(:); P=fft(S.ped); P=P(1:nf); np=size(S.d1,2);
        xp=p.isv(:)'/p.n; inr=xp>=XLO & xp<=XHI;
        Am=zeros(nf,np);
        for i=1:np, D=fft(S.d1(:,i)); Am(:,i)=abs(D(1:nf)./max(abs(P),eps)); end
        [~,jf]=min(abs(f-2)); e2=Am(jf,:); e2(~isfinite(e2))=0; e2(~inr)=0;
        fprintf('  gam=%d : peak@2k = %8.2f dB   (linear %.3e)\n', og, 20*log10(max(e2)), max(e2));
    end
    % BF sequence over the 7 standard places to see the fold directly
    pr=pa; pr.isv=[1136 1005 840 655 466 273 80];
    evalc('S=tdm26(0,pr,0,0);');
    nf=numel(S.f); f=S.f(:); P=fft(S.ped); P=P(1:nf); bf=nan(1,7);
    for i=1:7
        D=fft(S.d1(:,i)); H=D(1:nf)./max(abs(P),eps); a=abs(H).*f; a(~(f>0.15&f<18))=0;
        [pk,ip]=max(a); if(pk>0), bf(i)=f(ip); end
    end
    fprintf('  BF over 7 places (apex->base): %s kHz\n', num2str(bf,'%7.2f'));
end
disp('C3AMP_DIAG_DONE');
