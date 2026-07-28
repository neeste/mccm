% LEVEL-DEPENDENCE IN THE 4-CHAMBER, ON THE CORRECT AMPLITUDE KNOB.
%
% SN's prediction: as the amplitude of the OHC force decreases, the LATENCY and
% the BEST PLACE should both shift (tip collapses, peak moves basally, latency
% shortens) -- the mechanism behind the level-dependent frequency-place map.
%
% *** ohcgain, NOT gam, is the 4-chamber's OHC-amplitude knob. *** In the m>=4
% branch k_act = gh*k3 - gam*k4, so gam scales only the k4/r4 sub-term and the
% gh*k3 active coupling SURVIVES at gam=0. It is fsp=pa.ohcgain that scales the
% WHOLE force: act = sgn*fsp*(k_act.*d2 + r_act.*v2). Four earlier sweeps
% targeted gam=0.5, which is the wrong quantity; this is the right one and has
% never been swept for divergence.
%
% Everything else stays NATIVE (gam=1, chsz(4)=0.05, default sgn=+1) -- the
% native operating point is known sound (contrast 13.5 dB, xbest 0.400, lat
% 2.07 ms, degen 0), so this measures the real model, not a patched one.
%
% Stability on maxRe(all) (maxRe_osc is blind to the f=0 mode). Interior scoring
% x/L in [0.05,0.85]. Pre-flight gate first, per the standing rule.

gv  = [1.00 0.85 0.70 0.50 0.30 0.15 0.00];
fp  = [1 2 4 8];
ISV = 1391:-10:11; XLO=0.05; XHI=0.85;

mkpa = @(g) setfield(modpar26(4),'ohcgain',g); %#ok<SFLD>

% ---- pre-flight: does ohcgain move the output at all? ----
fprintf('\n=== PRE-FLIGHT: ohcgain 1 vs 0 ===\n');
pk=nan(1,2); gg=[1 0];
for j=1:2
    t=local_tip(mkpa(gg(j)),ISV,XLO,XHI,fp);
    pk(j)=t.pk(2); fprintf('  ohcgain=%.2f  peak@2k = %8.2f dB   contrast %6.1f\n', gg(j), t.pk(2), t.chi);
end
if (all(isfinite(pk)) && abs(pk(1)-pk(2))<1e-9)
    fprintf('  *** ABORT: ohcgain does not move the output.\n'); disp('OHCGAIN_ABORTED'); return
end
fprintf('  -> ohcgain acts (delta %.2f dB). Proceeding.\n', abs(pk(1)-pk(2)));

T=cell(numel(gv),1);
for i=1:numel(gv), T{i}=local_tip(mkpa(gv(i)),ISV,XLO,XHI,fp); end

fprintf('\n=========== BEST PLACE x/L vs ohcgain (basal = smaller) ===========\n');
fprintf('  ohcgain |  1 kHz   2 kHz   4 kHz   8 kHz\n');
for i=1:numel(gv)
    fprintf('   %5.2f   | %6.3f  %6.3f  %6.3f  %6.3f\n', gv(i), T{i}.xb);
end
fprintf('\n=========== LATENCY (ms) vs ohcgain ===========\n');
fprintf('  ohcgain |  1 kHz   2 kHz   4 kHz   8 kHz\n');
for i=1:numel(gv)
    fprintf('   %5.2f   | %6.2f  %6.2f  %6.2f  %6.2f\n', gv(i), T{i}.lt);
end
fprintf('\n=========== TIP / GAIN / STABILITY vs ohcgain ===========\n');
fprintf('  ohcgain  finite  contrast(dB)  peak@2k(dB)  gain(dB re 0)  deg\n');
p0=T{end}.pk(2);
for i=1:numel(gv)
    fprintf('   %5.2f   %-6s  %9.1f   %10.1f   %+11.1f    %d\n', ...
        gv(i), T{i}.fin, T{i}.chi, T{i}.pk(2), T{i}.pk(2)-p0, T{i}.deg);
end

fprintf('\n=========== coupeig stability vs ohcgain (native gam=1) ===========\n');
fprintf('  ohcgain   maxRe(all)    maxRe_osc    verdict\n');
for i=1:numel(gv)
    pa=mkpa(gv(i));
    try
        evalc('E=tdm26(''coupeig'',struct(''pa'',pa));');
        vv='STABLE (at floor)'; if (E.maxRe>1), vv='STATIC-UNSTABLE'; end
        if (E.maxRe_osc>0), vv='OSC-UNSTABLE'; end
        fprintf('   %5.2f   %+11.1f  %+10.1f    %s\n', gv(i), E.maxRe, E.maxRe_osc, vv);
    catch e
        fprintf('   %5.2f   coupeig failed: %s\n', gv(i), e.message);
    end
end
fprintf(['\nSN PREDICTION CONFIRMED if, as ohcgain falls: gain and contrast DROP,\n' ...
         'best place moves BASALLY (x/L smaller), latency SHORTENS -- all together,\n' ...
         'with degen 0 and maxRe(all) at the floor throughout.\n']);
disp('OHCGAIN_SWEEP_DONE');

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
