% m5 SWEEP AT chsz(4)=0.025  (corrected design)
%
% The first attempt swept m5 at chsz(4)=0 and returned bit-for-bit identical
% rows. m5o was NOT inert -- cp.m5 and mu3 changed exactly as set. The cause was
% an EXACT ALGEBRAIC CANCELLATION: with chsz(4)=0, L4_p=L4_c=0, so the CL row
% (tdm26.m:422) is a2(k,14)=-mu3, a2(k,16)=+mu3, a1/a3(k,16)=0, against RHS
% qx(j4)=s3.*mu3 (tdm26.m:603):
%       -mu3*p14 + mu3*p16 = mu3*s3   ->   p16 - p14 = s3
% mu3 divides out. Zeroing the CL fluid inertance makes the CL chamber
% independent of m5, so chsz(4) and m5 are NOT independent levers.
%
% Here chsz(4)=0.025 (half native 0.05): the fluid arm is reduced but nonzero,
% so mu3 retains purchase. Reference from cl_coupling.m: chsz(4)=0.025, m5 x1
% gave maxRe(all)@gam=0.5 = +15685.8.
%
% PRE-FLIGHT GATE (lesson from the wasted sweep): confirm m5 moves maxRe ITSELF
% -- not merely the intermediate mu3, which moved correctly last time and still
% told us nothing -- before spending the full sweep. 2 runs instead of 7 if flat.

mult = [0.5 1 2 4 8 16 32];
m5b  = 0.0360276; k5b = 3.05084e8; CH4 = 0.025;
ISV  = 1391:-10:11; XLO=0.05; XHI=0.85;

mk = @(g,mm) setfield(setfield(setfield(modpar26(4),'gam',g), ...
                     'chsz',[0.95 0.05 1.0 CH4]),'m5o',m5b*mm); %#ok<SFLD>

fprintf('\n=== PRE-FLIGHT: does m5 move maxRe at chsz(4)=%.3f? ===\n', CH4);
pf=nan(1,2); pv=[1 8];
for j=1:2
    pa=mk(0.5,pv(j));
    try
        evalc('E=tdm26(''coupeig'',struct(''pa'',pa));');
        pf(j)=E.maxRe;
        fprintf('  m5 x%-3g  maxRe(all) = %+11.1f   maxRe_osc = %+9.1f\n', pv(j), E.maxRe, E.maxRe_osc);
    catch e
        fprintf('  m5 x%-3g  FAILED: %s\n', pv(j), e.message);
    end
end
if (all(isfinite(pf)) && abs(pf(1)-pf(2)) < 1e-6)
    fprintf('\n  *** ABORT: maxRe IDENTICAL across m5 x1 vs x8 at chsz(4)=%.3f.\n', CH4);
    fprintf('  m5 still has no influence on the output -- do NOT run the full sweep.\n');
    disp('CL_M5B_ABORTED'); return
end
fprintf('  -> m5 DOES move maxRe (delta = %.1f). Proceeding.\n', abs(pf(1)-pf(2)));

fprintf('\n=========== m5 sweep at chsz(4)=%.3f, m=4 ===========\n', CH4);
fprintf('  m5(x)  mu3rel  fOC(kHz) | gam=1.0 tip (guardrail)      | gam=0.5\n');
fprintf('                          | fin    contrast xbest  lat deg| fin    deg\n');
fprintf('%s\n', repmat('-',1,86));
for im=1:numel(mult)
    mm=mult(im); foc=sqrt(k5b/(m5b*mm))/2/pi/1000;
    t1=local_tip(mk(1.0,mm),ISV,XLO,XHI); t5=local_tip(mk(0.5,mm),ISV,XLO,XHI);
    fprintf('  %5.1f  %6.3f  %7.2f | %-6s %7.1f %6.3f %5.2f %d | %-6s %d\n', ...
        mm, 1/mm, foc, t1.fin, t1.chi, t1.xb, t1.lt, t1.deg, t5.fin, t5.deg);
end

fprintf('\n=========== coupeig maxRe(all) at gam=0.5, chsz(4)=%.3f ===========\n', CH4);
fprintf('  m5(x)  mu3rel    maxRe(all)    maxRe_osc    verdict\n');
for im=1:numel(mult)
    mm=mult(im); pa=mk(0.5,mm);
    try
        evalc('E=tdm26(''coupeig'',struct(''pa'',pa));');
        vv='STABLE (at floor)'; if (E.maxRe>1), vv='STATIC-UNSTABLE'; end
        if (E.maxRe_osc>0), vv='OSC-UNSTABLE'; end
        fprintf('  %5.1f  %6.3f  %+11.1f  %+10.1f    %s\n', mm, 1/mm, E.maxRe, E.maxRe_osc, vv);
    catch e
        fprintf('  %5.1f  %6.3f  coupeig failed: %s\n', mm, 1/mm, e.message);
    end
end
fprintf('\nReference: chsz(4)=0.025, m5 x1 gave maxRe(all)=+15685.8 (cl_coupling.m).\n');
fprintf('Target: SMALLEST m5 reaching the +0.0 floor with gam=1 deg=0 and contrast intact.\n');
disp('CL_M5B_DONE');

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
