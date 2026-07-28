% CL<->SS COUPLING LEVER: does reducing the CL fluid return path break the
% basal static feedback loop (shear -> active force -> OC-height -> CL<->SS ->
% shear) that the eigenvector decomposition localized to OC-height at the base?
%
% chsz(4) is the CL chamber's longitudinal FLUID inertance (tdm26.m:412 L4;
% enters the CL row a2(k,16)=L4_p+L4_c+mu3 and the longitudinal a1/a3(k,16)).
% It is the RETURN arm of the loop. (Distinct from mu3=m1/m5, the direct SS<->CL
% MECHANICAL mass coupling; noted, not swept here.)
%
% NOTE tdm26.m:353 renormalizes pa.chsz to sum=2, so setting chsz(4)=v rescales
% the others. Native chsz=[0.95 0.05 1.0 0.05] (sum 2.05); effective chsz(4) =
% v*2/(2.0+v). Both raw and effective are reported.
%
% PASS = reducing chsz(4) drives maxRe(all) at gam=0.5 toward the +0.0 floor
% (loop broken) WITHOUT flattening the native gam=1 tip. Stability read on
% maxRe(all) (not maxRe_osc, blind to this f=0 mode). Interior scoring for tips.

v4   = [0 0.01 0.025 0.05 0.10 0.20];   % raw chsz(4); 0.05 is native
ISV  = 1391:-10:11; XLO=0.05; XHI=0.85;

clicktip = @(pa) local_tip(pa, ISV, XLO, XHI);

fprintf('\n=========== CL fluid coupling chsz(4) sweep, m=4 ===========\n');
fprintf('  chsz4  eff4  | gam=1.0 native tip           | gam=0.5 collapse\n');
fprintf('  (raw)       | finite contrast xbest  lat deg| finite contrast xbest  lat deg\n');
fprintf('%s\n', repmat('-',1,92));
for iv=1:numel(v4)
    v=v4(iv); eff=v*2/(2.0+v);
    pa1=modpar26(4); pa1.gam=1.0; pa1.chsz(4)=v;
    pa5=modpar26(4); pa5.gam=0.5; pa5.chsz(4)=v;
    t1=clicktip(pa1); t5=clicktip(pa5);
    fprintf('  %5.3f %5.3f| %-6s %7.1f %6.3f %5.2f %d | %-6s %7.1f %6.3f %5.2f %d\n', ...
        v, eff, t1.fin, t1.chi, t1.xb, t1.lt, t1.deg, ...
                t5.fin, t5.chi, t5.xb, t5.lt, t5.deg);
end

fprintf('\n=========== coupeig maxRe(all) at gam=0.5 vs chsz(4) ===========\n');
fprintf('  chsz4(raw)  eff4    maxRe(all)    maxRe_osc    verdict\n');
for iv=1:numel(v4)
    v=v4(iv); eff=v*2/(2.0+v);
    pa=modpar26(4); pa.gam=0.5; pa.chsz(4)=v;
    try
        evalc('E=tdm26(''coupeig'',struct(''pa'',pa));');
        vv='stable-ish'; if (E.maxRe>1), vv='STATIC-UNSTABLE'; end
        if (E.maxRe_osc>0), vv='OSC-UNSTABLE'; end
        fprintf('  %6.3f     %5.3f  %+11.1f  %+10.1f    %s\n', v, eff, E.maxRe, E.maxRe_osc, vv);
    catch e
        fprintf('  %6.3f     %5.3f  coupeig failed: %s\n', v, eff, e.message);
    end
end
fprintf('\nPASS: maxRe(all)@0.5 falls toward +0.0 as chsz(4) drops, gam=1 tip intact.\n');
disp('CL_COUPLING_DONE');

function t=local_tip(pa, ISV, XLO, XHI)
t.fin='yes'; t.chi=NaN; t.xb=NaN; t.lt=NaN; t.deg=1;
pa.isv=ISV;
try, evalc('S=tdm26(0,pa,0,0);'); catch, t.fin='THREW'; return; end
if (any(~isfinite(S.d1(:)))), t.fin='DIVERG'; return; end
nf=numel(S.f); f=S.f(:); P=fft(S.ped); P=P(1:nf); np=size(S.d1,2);
xp=pa.isv(:)'/pa.n; inr=xp>=XLO & xp<=XHI;
H=zeros(nf,np); for i=1:np, D=fft(S.d1(:,i)); H(:,i)=D(1:nf)./max(abs(P),eps); end
Am=abs(H); co=nan(1,np);
for i=find(inr)
    a=Am(:,i).*f; a(~(f>0.15&f<18))=0; [pk,ip]=max(a); if(pk<=0),continue;end
    if(f(ip)<0.5||f(ip)>16.5),continue;end
    [~,it]=min(abs(f-f(ip)/4)); co(i)=20*log10(pk/max(a(it),eps));
end
ii=find(isfinite(co)); if(~isempty(ii)), t.chi=co(ii(end)); end
w=2*pi*f*1000; [~,ifq]=min(abs(f-2)); ex=Am(ifq,:); ex(~isfinite(ex))=0; ex(~inr)=0;
[~,ib]=max(ex);
if(ex(ib)>0)
    t.xb=xp(ib); ph=unwrap(angle(H(:,ib))); gd=-gradient(ph,w);
    t.lt=median(gd(max(1,ifq-3):min(nf,ifq+3)))*1000;
    XB=nan(1,4); LT=XB; fp=[1 2 4 8];
    for k=1:4
        [~,jf]=min(abs(f-fp(k))); e2=Am(jf,:); e2(~isfinite(e2))=0; e2(~inr)=0;
        [~,jb]=max(e2); if(e2(jb)<=0),continue;end
        XB(k)=xp(jb); p2=unwrap(angle(H(:,jb))); g2=-gradient(p2,w);
        LT(k)=median(g2(max(1,jf-3):min(nf,jf+3)))*1000;
    end
    t.deg=(std(XB,0,'omitnan')<1e-3)||(median(abs(LT),'omitnan')<0.30);
end
end
