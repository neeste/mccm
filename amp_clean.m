% CLEAN amplifier re-measurement on the decompressing configs.
% mass_cost.m's gain column was untrustworthy (native gain@2k read +0.7 dB vs
% the dedicated ohcgain sweep's +2.4 dB). This uses the COMPLETE, validated
% local_tip from ohcgain_sweep.m and a native CONTROL that must reproduce +2.4.
%
% Question: does decompressing the map via heavy mass preserve the OHC AMPLIFIER
% (gain re passive) even though it flattens the passive CONTRAST (14.7 -> 2 dB)?
ISV=1391:-10:11; XLO=0.05; XHI=0.85; fp=[1 2 4 8];
b2=modpar26(4).m2o; b5=modpar26(4).m5o;
cfg={ {'native', 1, 1}, {'m2o x32', 32, 1}, {'m5o x16', 1, 16}, {'m5o x32', 1, 32} };

fprintf('\n  config    |  peak@2k(dB): og=1 / og=0 | gain(1vs0)@2k | contrast(og=1)\n');
fprintf('%s\n', repmat('-',1,72));
for c=1:numel(cfg)
    nm=cfg{c}{1}; m2m=cfg{c}{2}; m5m=cfg{c}{3};
    pa1=modpar26(4); pa1.m2o=b2*m2m; pa1.m5o=b5*m5m; pa1.ohcgain=1;
    pa0=modpar26(4); pa0.m2o=b2*m2m; pa0.m5o=b5*m5m; pa0.ohcgain=0;
    t1=local_tip(pa1,ISV,XLO,XHI,fp);
    t0=local_tip(pa0,ISV,XLO,XHI,fp);
    g=t1.pk(2)-t0.pk(2);
    fprintf('  %-8s  |   %8.2f / %8.2f    |    %+6.2f     |   %6.1f\n', ...
            nm, t1.pk(2), t0.pk(2), g, t1.chi);
end
fprintf('%s\n', repmat('-',1,72));
fprintf('CONTROL: native gain(1vs0)@2k must be ~+2.4 dB (matches ohcgain_sweep).\n');
fprintf('If it is, read the rest: amplifier SURVIVES decompression if the heavy-\n');
fprintf('mass configs keep gain near native despite their collapsed contrast.\n');
disp('AMP_CLEAN_DONE');

function t=local_tip(pa, ISV, XLO, XHI, fp)
t.fin='yes'; t.chi=NaN; t.deg=1; t.xb=nan(1,numel(fp)); t.lt=t.xb; t.pk=t.xb;
pa.isv=ISV;
try, evalc('S=tdm26(0,pa,0,0);'); catch, t.fin='THREW'; return; end
if (any(~isfinite(S.d1(:)))), t.fin='DIVERG'; return; end
nf=numel(S.f); f=S.f(:); P=fft(S.ped); P=P(1:nf); np=size(S.d1,2);
xp=pa.isv(:)'/pa.n; inr=xp>=XLO & xp<=XHI;
H=zeros(nf,np); for i=1:np, D=fft(S.d1(:,i)); H(:,i)=D(1:nf)./max(abs(P),eps); end
Am=abs(H); co=nan(1,np);
for i=find(inr)
    a=Am(:,i).*f; a(~(f>0.15&f<18))=0; [q,ip]=max(a); if(q<=0),continue;end
    if(f(ip)<0.5||f(ip)>16.5),continue;end
    [~,it]=min(abs(f-f(ip)/4)); co(i)=20*log10(q/max(a(it),eps));
end
ii=find(isfinite(co)); if(~isempty(ii)), t.chi=co(ii(end)); end
w=2*pi*f*1000;
for k=1:numel(fp)
    [~,jf]=min(abs(f-fp(k))); e2=Am(jf,:); e2(~isfinite(e2))=0; e2(~inr)=0;
    [~,jb]=max(e2); if(e2(jb)<=0),continue;end
    t.xb(k)=xp(jb); t.pk(k)=20*log10(e2(jb));
    p2=unwrap(angle(H(:,jb))); g2=-gradient(p2,w);
    t.lt(k)=median(g2(max(1,jf-3):min(nf,jf+3)))*1000;
end
t.deg=(std(t.xb,0,'omitnan')<1e-3)||(median(abs(t.lt),'omitnan')<0.30);
end
