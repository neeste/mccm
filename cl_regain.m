% RE-TEST OF THE CL AMPLIFICATION LIMIT, WITH THE -gam*r4 FIX.
%
% The standing conclusion -- "CL cannot amplify: 1.40x max stable, physically
% unstable above gain ~3, coupeig +470.9" -- rests on two things that are now
% both in doubt:
%   (1) the r4 bug: the m>=4 r_act was missing -gam*cp.r4, so the OHC force was
%       a pure negative stiffness and could not inject energy at ANY gain. The
%       "instability above gain ~3" was STATIC DIVERGENCE from net negative
%       stiffness, not a physical amplification ceiling.
%   (2) the metric: cl_gain.m / cl_gain2.m both scored amplification with
%       max(S.wnr). WNR magnitude SATURATES -- it must not be used to infer
%       amplifier behavior. The 1.40x figure inherits that flaw regardless of (1).
%
% So this re-tests with three independent metrics:
%   A  linearized click gain (dB re ohcgain=0) -- reliable amplitude measure,
%      scored on the interior only (x/L in [0.05,0.85]); includes r4act=0 control
%   B  ohcP energy diagnostic -- the physical arbiter: >0 = force opposes damping
%      and injects energy. Requires the wnr1 path (dgn is built there).
%   C  coupeig -- where the true stability ceiling now sits.

GV = [0 0.3 1.0 2.0 3.0];
SG = [+1 -1];
ISV = 1391:-10:11; XLO=0.05; XHI=0.85;

% ================= A: linearized click gain (fast, reliable) =================
fprintf('\n===== A: linearized click gain, m=4 (interior x/L %.2f-%.2f) =====\n',XLO,XHI);
fprintf(' r4act  sgn  ohcgain   peak@2k(dB)  gain(dB)   contrast(dB)  lat@2k(ms)\n');
for ab = [0 1]
for si = 1:numel(SG)
    pk0 = NaN;
    for gi = 1:numel(GV)
        pa = modpar26(4); pa.isv = ISV; pa.r4act = ab;
        pa.ohcsgn = SG(si); pa.ohcgain = GV(gi);
        try
            evalc('S = tdm26(0,pa,0,0);');
        catch e
            fprintf(' %5d  %+3d   %5.2f    FAILED: %s\n', ab,SG(si),GV(gi),e.message); continue
        end
        nf=numel(S.f); f=S.f(:); P=fft(S.ped); P=P(1:nf);
        np=size(S.d1,2); xp=pa.isv(:)'/pa.n; inr = xp>=XLO & xp<=XHI;
        H=zeros(nf,np);
        for i=1:np, D=fft(S.d1(:,i)); H(:,i)=D(1:nf)./max(abs(P),eps); end
        Am=abs(H);
        co=nan(1,np);
        for i=find(inr)
            a=Am(:,i).*f; a(~(f>0.15&f<18))=0;
            [pk,ip]=max(a); if (pk<=0), continue; end
            if (f(ip)<0.5||f(ip)>16.5), continue; end
            [~,it]=min(abs(f-f(ip)/4)); co(i)=20*log10(pk/max(a(it),eps));
        end
        iv=find(isfinite(co)); chi=NaN; if (~isempty(iv)), chi=co(iv(end)); end
        [~,ifq]=min(abs(f-2)); ex=Am(ifq,:); ex(~isfinite(ex))=0; ex(~inr)=0;
        [~,ib]=max(ex); pkd=NaN; lt=NaN;
        if (ex(ib)>0)
            pkd=20*log10(ex(ib));
            w=2*pi*f*1000; ph=unwrap(angle(H(:,ib))); gd=-gradient(ph,w);
            lt=median(gd(max(1,ifq-3):min(nf,ifq+3)))*1000;
        end
        if (gi==1), pk0=pkd; end
        fprintf(' %5d  %+3d   %5.2f    %9.1f   %+7.1f      %7.1f     %7.2f\n', ...
                ab, SG(si), GV(gi), pkd, pkd-pk0, chi, lt);
    end
end
end

% ================= B: energy diagnostic (the physical arbiter) ===============
fprintf('\n===== B: ohcP energy diagnostic, r4act=1, wnr1 @ 2 kHz, 20 dB =====\n');
fprintf('  sgn  ohcgain   lat(ms)   max|WNR|        ohcP          ohcW      nonfin\n');
for si = 1:numel(SG)
    for gi = [1 3 5]                       % ohcgain = 0, 1.0, 3.0
        p.fr=2; p.lv=20; p.pa=modpar26(4); p.pa.hbmode='bm';
        p.pa.ohcsgn=SG(si); p.pa.ohcgain=GV(gi); p.pa.r4act=1;
        try
            evalc('S=tdm26(''wnr1'',p,0,0);'); d=S.dgn;
            fprintf('  %+3d   %5.2f   %7.2f  %11.3e %+12.4e %+12.4e %6d\n', ...
                    SG(si),GV(gi),S.tpk,max(S.wnr),d.ohcP,d.ohcW,d.ohcNaN);
        catch e
            fprintf('  %+3d   %5.2f   FAILED: %s\n', SG(si),GV(gi),e.message);
        end
    end
end
fprintf('  (ohcgain=0 is the validity control: ohcP must be ~0)\n');

% ================= C: where the stability ceiling actually is ===============
fprintf('\n===== C: coupeig stability ceiling, r4act=1 =====\n');
fprintf('  sgn  ohcgain   maxRe_osc   verdict\n');
for si = 1:numel(SG)
    for g = [1.0 2.0 3.0 5.0]
        pa=modpar26(4); pa.ohcsgn=SG(si); pa.ohcgain=g; pa.r4act=1;
        try
            evalc('E=tdm26(''coupeig'',struct(''pa'',pa));');
            v='sub-critical (stable)'; if (E.maxRe_osc>=0), v='UNSTABLE'; end
            fprintf('  %+3d   %5.2f   %+9.1f   %s\n', SG(si), g, E.maxRe_osc, v);
        catch e
            fprintf('  %+3d   %5.2f   coupeig failed: %s\n', SG(si), g, e.message);
        end
    end
end
disp('CL_REGAIN_DONE');
