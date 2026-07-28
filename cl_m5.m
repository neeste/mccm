% BREAK THE REMAINING LOOP: sweep m5 at chsz(4)=0.
%
% The CL fluid return arm is already cut (chsz(4)=0 dropped maxRe(all)@gam=0.5
% from +16862 to +3016, an 82% cut, with the gam=1 tip intact). The residual
% +3016 comes through the MECHANICAL return arm mu3 = cp.m1./cp.m5 (tdm26.m:407),
% the direct SS<->CL mass coupling that the unstable eigenvector implicated
% (OC-height 69% coupled to shear 26%).
%
% mu3 = m1/m5, so RAISING m5o CUTS the mechanical return coupling. cp.m5 =
% m5o*exp(m5e*x+m5q*q) (tdm26.m:477); modpar26c4 m5o = 0.0360276.
%
% WATCH (m5 is not a free lever): larger m5 also
%   (a) lowers the OC-height resonance ~sqrt(k5/m5) -- base is
%       sqrt(3.05e8/0.036)/2pi ~ 14.6 kHz, so x16 puts it ~3.7 kHz, IN BAND;
%   (b) reduces OC-height participation via a(i3) = a(i1) + s3./cp.m5 (line 630).
% So it may cut the TIP along with the loop. Target = the SMALLEST m5 that
% drives maxRe(all) below the +0.0 floor while keeping the gam=1 tip (degen 0).
% mult=0.5 is the direction check: it should make the instability WORSE.
%
% Stability read on maxRe(all). NOTE from the chsz sweep: a "finite" click is
% NOT proof of stability (the basal static mode can be weakly excited), so the
% click columns are the TIP guardrail, not the stability verdict.

mult = [0.5 1 2 4 8 16 32];
m5b  = 0.0360276; k5b = 3.05084e8;
ISV  = 1391:-10:11; XLO=0.05; XHI=0.85;

fprintf('\n=========== m5 sweep at chsz(4)=0, m=4 ===========\n');
fprintf('  m5(x)  mu3rel  fOC(kHz) | gam=1.0 tip (guardrail)     | gam=0.5\n');
fprintf('                          | fin    contrast xbest  lat deg| fin    deg\n');
fprintf('%s\n', repmat('-',1,86));
for im=1:numel(mult)
    mm=mult(im); foc=sqrt(k5b/(m5b*mm))/2/pi/1000;
    pa1=modpar26(4); pa1.gam=1.0; pa1.chsz(4)=0; pa1.m5o=m5b*mm;
    pa5=modpar26(4); pa5.gam=0.5; pa5.chsz(4)=0; pa5.m5o=m5b*mm;
    t1=local_tip(pa1,ISV,XLO,XHI); t5=local_tip(pa5,ISV,XLO,XHI);
    fprintf('  %5.1f  %6.3f  %7.2f | %-6s %7.1f %6.3f %5.2f %d | %-6s %d\n', ...
        mm, 1/mm, foc, t1.fin, t1.chi, t1.xb, t1.lt, t1.deg, t5.fin, t5.deg);
end

fprintf('\n=========== coupeig maxRe(all) at gam=0.5, chsz(4)=0 ===========\n');
fprintf('  m5(x)  mu3rel    maxRe(all)    maxRe_osc    verdict\n');
for im=1:numel(mult)
    mm=mult(im);
    pa=modpar26(4); pa.gam=0.5; pa.chsz(4)=0; pa.m5o=m5b*mm;
    try
        evalc('E=tdm26(''coupeig'',struct(''pa'',pa));');
        vv='STABLE (at floor)'; if (E.maxRe>1), vv='STATIC-UNSTABLE'; end
        if (E.maxRe_osc>0), vv='OSC-UNSTABLE'; end
        fprintf('  %5.1f  %6.3f  %+11.1f  %+10.1f    %s\n', mm, 1/mm, E.maxRe, E.maxRe_osc, vv);
    catch e
        fprintf('  %5.1f  %6.3f  coupeig failed: %s\n', mm, 1/mm, e.message);
    end
end
fprintf('\nReference: chsz(4)=0, m5 x1 gave maxRe(all)=+3016 (from cl_coupling.m).\n');
fprintf('Looking for the SMALLEST m5 with maxRe(all) at the +0.0 floor AND gam=1 deg=0.\n');
disp('CL_M5_DONE');

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
