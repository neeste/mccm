% DOES A STIFFER OC-HEIGHT COORDINATE REMOVE THE m=4 gam=0.5 COLLAPSE?
% (Distinct from the earlier k5_sweep.m, which sought the ohcP energy-injection
%  phase via ohcgain/ohcsgn. Here gam is the OHC amplitude and the question is
%  STATIC STABILITY of the passive OC-height resonance.)
%
% SN's design principle: the OHC force pushes BETWEEN BM and RL, an internal
% force pair on the OC-height coordinate d3 (tdm26.m:727 s3 = -(k5*dc + r5*vc +
% act)). An internal force pair does no net positional constraint, so d3 needs
% its OWN passive restoring stiffness cp.k5 to be a bounded resonator; too soft
% a k5 leaves a positive real (STATIC) eigenvalue -> divergence. cp.k5 = k5o *
% exp(k5e*x + k5q*q) (tdm26.m:475), so k5o scales the restoring force uniformly.
%
% The provisional modpar26c4 value is k5o = 3.05e8. Earlier, m=4 collapsed at
% gam=0.5 (-64.9 dB pinned at the basal edge, zero latency, degenerate). If that
% collapse is under-restoration, RAISING k5o should convert it to a clean
% interior tip. If k5o does nothing, the collapse is something else.
%
% Stability is read as TIME-MARCH FINITENESS (the ground truth that maxRe_osc
% missed) plus coupeig maxRe(ALL modes) -- NOT maxRe_osc, which cannot see the
% static f=0 divergence responsible here. Interior scoring x/L in [0.05,0.85].

mult = [0.5 1 2 4 8];
base = 3.05084e8;
gv   = [0.3 0.5 0.7 1.0];
ISV  = 1391:-10:11; XLO=0.05; XHI=0.85;

% ---------- PART 1: click sweep (fast; divergence + tip quality) ----------
fprintf('\n=========== m=4 click sweep: k5o multiplier x gam ===========\n');
for mi = 1:numel(mult)
    fprintf('\n--- k5o = %.3g  (x%.2g) ---\n', base*mult(mi), mult(mi));
    fprintf('  gam   finite    contrast(dB) peak@2k(dB)  xbest@2k  lat@2k(ms)  degen\n');
    for gi = 1:numel(gv)
        pa = modpar26(4); pa.gam = gv(gi); pa.isv = ISV;
        pa.k5o = base*mult(mi);
        fin='yes'; d1bad=false;
        try, evalc('S = tdm26(0,pa,0,0);'); catch e
            fprintf('  %4.2f   THREW: %s\n', gv(gi), e.message); continue; end
        if (any(~isfinite(S.d1(:)))), fin='DIVERGES'; d1bad=true; end
        nf=numel(S.f); f=S.f(:); P=fft(S.ped); P=P(1:nf);
        np=size(S.d1,2); xp=pa.isv(:)'/pa.n; inr = xp>=XLO & xp<=XHI;
        chi=NaN; pkd=NaN; xb=NaN; lt=NaN; deg=1;
        if (~d1bad)
            H=zeros(nf,np);
            for i=1:np, D=fft(S.d1(:,i)); H(:,i)=D(1:nf)./max(abs(P),eps); end
            Am=abs(H); co=nan(1,np);
            for i=find(inr)
                a=Am(:,i).*f; a(~(f>0.15&f<18))=0;
                [pk,ip]=max(a); if (pk<=0), continue; end
                if (f(ip)<0.5||f(ip)>16.5), continue; end
                [~,it]=min(abs(f-f(ip)/4)); co(i)=20*log10(pk/max(a(it),eps));
            end
            iv=find(isfinite(co)); if (~isempty(iv)), chi=co(iv(end)); end
            [~,ifq]=min(abs(f-2)); ex=Am(ifq,:); ex(~isfinite(ex))=0; ex(~inr)=0;
            [~,ib]=max(ex);
            if (ex(ib)>0)
                pkd=20*log10(ex(ib)); xb=xp(ib);
                w=2*pi*f*1000; ph=unwrap(angle(H(:,ib))); gd=-gradient(ph,w);
                lt=median(gd(max(1,ifq-3):min(nf,ifq+3)))*1000;
                XBv=nan(1,4); LTv=XBv; fp=[1 2 4 8];
                for k=1:4
                    [~,jf]=min(abs(f-fp(k))); e2=Am(jf,:); e2(~isfinite(e2))=0; e2(~inr)=0;
                    [~,jb]=max(e2); if (e2(jb)<=0), continue; end
                    XBv(k)=xp(jb); p2=unwrap(angle(H(:,jb))); g2=-gradient(p2,w);
                    LTv(k)=median(g2(max(1,jf-3):min(nf,jf+3)))*1000;
                end
                deg = (std(XBv,0,'omitnan')<1e-3) || (median(abs(LTv),'omitnan')<0.30);
            end
        end
        fprintf('  %4.2f   %-8s  %8.1f    %8.1f   %7.4f   %8.2f    %d\n', ...
                gv(gi), fin, chi, pkd, xb, lt, deg);
    end
end

% ---------- PART 2: coupeig maxRe(ALL) at gam=0.5 across k5o ----------
fprintf('\n=========== m=4 coupeig at gam=0.50 (static stability) ===========\n');
fprintf('  (maxRe(all) is the true indicator; maxRe_osc is blind to f=0 divergence)\n');
fprintf('  k5o(x)    maxRe(all)     maxRe_osc     verdict\n');
for mi = 1:numel(mult)
    pa=modpar26(4); pa.gam=0.5; pa.k5o=base*mult(mi);
    try
        evalc('E=tdm26(''coupeig'',struct(''pa'',pa));');
        v='stable'; if (E.maxRe > 1), v='STATIC-UNSTABLE'; end
        if (E.maxRe_osc > 0), v='OSC-UNSTABLE'; end
        fprintf('  %5.2f   %+11.1f  %+11.1f    %s\n', mult(mi), E.maxRe, E.maxRe_osc, v);
    catch e
        fprintf('  %5.2f   coupeig failed: %s\n', mult(mi), e.message);
    end
end
fprintf('\nSN prediction: raising k5o removes the gam=0.5 collapse (finite, degen 0,\n');
fprintf('interior tip) and drives maxRe(all) from + toward the +0.0 stable floor.\n');
disp('K5_COLLAPSE_DONE');
