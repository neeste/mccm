% Does decompressing the map via heavy mass COST the amplifier and/or stability?
% Hypothesis: heavy m2/m5 decouples the extra DOF -> wider map, but that DOF is
% where the graded OHC amplifier lives, so the amp may collapse. Check the best
% decompressing configs on THREE axes: range (recompute), stability maxRe(all),
% and amplifier = gain re passive (ohcgain 1 vs 0) + contrast at gam=1.
ISV=[1136 1005 840 655 466 273 80]; b2=modpar26(4).m2o; b5=modpar26(4).m5o;
cfg={ {'native',        1,1}, {'m2o x32 (best rng)',32,1}, ...
      {'m5o x16',        1,16},{'m5o x32',           1,32} };
fprintf('\n  config              range  mono | maxRe(all)  osc | contrast  gain(1vs0 dB)\n');
for c=1:numel(cfg)
    nm=cfg{c}{1}; pa=modpar26(4); pa.m2o=b2*cfg{c}{2}; pa.m5o=b5*cfg{c}{3};
    % range + contrast at ohcgain=1
    [rng,mono]=rng_tip(pa,ISV,1); [~,~,chi1,pk1]=rng_tip(pa,ISV,1);
    [~,~,~,pk0]=rng_tip(pa,ISV,0);          % passive
    gaindb=pk1-pk0;
    % stability
    mr=NaN; mo=NaN;
    try, evalc('E=tdm26(''coupeig'',struct(''pa'',pa));'); mr=E.maxRe; mo=E.maxRe_osc; catch, end
    fprintf('  %-18s  %5.2f  %-4s | %+9.1f %+7.1f | %6.1f    %+6.1f\n', ...
            nm, rng, mono, mr, mo, chi1, gaindb);
end
fprintf('\nTENSION CONFIRMED if range up but gain(1vs0) collapses toward 0 -- the\n');
fprintf('extra DOF was frozen out, taking the amplifier with it.\n');
disp('MASS_COST_DONE');

function [rng,mono,chi,pk2]=rng_tip(pa,ISV,ohcg)
rng=NaN;mono='?';chi=NaN;pk2=NaN; pa.isv=ISV; pa.ohcgain=ohcg;
try, evalc('S=tdm26(0,pa,0,0);'); catch, return; end
if(any(~isfinite(S.d1(:)))),mono='DIVG';return;end
nf=numel(S.f);f=S.f(:);P=fft(S.ped);P=P(1:nf);bf=nan(1,numel(ISV));co=bf;
for i=1:numel(ISV)
    D=fft(S.d1(:,i));H=D(1:nf)./max(abs(P),eps);a=abs(H).*f;a(~(f>0.15&f<18))=0;
    [pk,ip]=max(a);if(pk>0),bf(i)=f(ip);[~,it]=min(abs(f-f(ip)/4));co(i)=20*log10(pk/max(a(it),eps));end
end
v=bf(isfinite(bf)&bf>0); if(numel(v)>=2),rng=log2(max(v)/min(v));end
bx=bf(isfinite(bf)); if(numel(bx)>=3),if(all(diff(bx)>0)),mono='ok';else,mono='FOLD';end;end
iv=find(isfinite(co)); if(~isempty(iv)),chi=co(iv(end));end
[~,jf]=min(abs(f-2)); e2=abs(H); % pk@2k on last-computed H is wrong; recompute properly below
% peak @2 kHz across places:
Am=zeros(nf,numel(ISV));
for i=1:numel(ISV), D=fft(S.d1(:,i)); Am(:,i)=abs(D(1:nf)./max(abs(P),eps)); end
xp=ISV/pa.n; inr=xp>=0.05&xp<=0.85; e2=Am(jf,:); e2(~inr)=0; pk2=20*log10(max(e2));
end
