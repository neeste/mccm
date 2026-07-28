ISV=1391:-10:11; XLO=0.05; XHI=0.85; fp=[1 2 4 8];
b2=modpar26(4).m2o;
L=load('parfit26_c4map_m2.mat'); paf=L.R.pa;
cfg={ {'m=2 target', modpar26(2)}, {'m=4 native', modpar26(4)}, ...
      {'m2o x32 (unfit)', setfield(modpar26(4),'m2o',b2*32)}, ...
      {'m2o x32 FITTED', paf} };
fprintf('\n  config           | range(oct) mono | contrast | gain(1vs0)@2k | maxRe\n');
fprintf('%s\n', repmat('-',1,74));
for c=1:numel(cfg)
    nm=cfg{c}{1}; pa=cfg{c}{2};
    pa1=pa; pa1.ohcgain=1; t1=local_tip(pa1,ISV,XLO,XHI,fp);
    pa0=pa; pa0.ohcgain=0; t0=local_tip(pa0,ISV,XLO,XHI,fp);
    % BF-map range over the 7 standard places
    rng=NaN; mono='?';
    pr=pa; pr.isv=[1136 1005 840 655 466 273 80];
    evalc('S=tdm26(0,pr,0,0);'); nf=numel(S.f); f=S.f(:); P=fft(S.ped); P=P(1:nf); bf=nan(1,7);
    for i=1:7, D=fft(S.d1(:,i)); H=D(1:nf)./max(abs(P),eps); a=abs(H).*f; a(~(f>0.15&f<18))=0; [pk,ip]=max(a); if(pk>0),bf(i)=f(ip);end; end
    v=bf(isfinite(bf)&bf>0); if(numel(v)>=2),rng=log2(max(v)/min(v));end
    bx=bf(isfinite(bf)); if(numel(bx)>=3),if(all(diff(bx)>0)),mono='ok';else,mono='FOLD';end;end
    mr=NaN; try, evalc('E=tdm26(''coupeig'',struct(''pa'',pa));'); mr=E.maxRe; catch, end
    fprintf('  %-16s |   %5.2f    %-4s |  %6.1f  |    %+6.2f     | %+.1f\n', ...
            nm, rng, mono, t1.chi, t1.pk(2)-t0.pk(2), mr);
end
fprintf('%s\n', repmat('-',1,74));
fprintf('  fitted profile: k1e=%.3f m1e=%.3f  chsz=[%s]  (m2o pinned %.3f)\n', ...
        paf.k1e, paf.m1e, num2str(paf.chsz,'%.2f '), paf.m2o);
disp('C4_FINAL_DONE');
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
