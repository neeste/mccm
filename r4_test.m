% A/B TEST OF THE MISSING gam*r4 ACTIVE-RESISTANCE TERM  (m = 3 and m = 4)
%
% fdm26 (the reference implementation) builds the active impedance as
%     z4 = k4/s + r4                                   (fdm26.m:960)
%     zg = gam*z4                                      (fdm26.m:852,869)
%     zk(1,1:2) = [z1, gh*zh - zg]                     (fdm26.m:890)
% so the row-1 coefficient of V2 is  gh*(k3/s+r3) - gam*(k4/s+r4), which in the
% time domain is
%     k_act = gh*k3 - gam*k4        <- tdm26 had this
%     r_act = gh*r3 - gam*r4        <- tdm26 DROPPED the -gam*r4
% cp.r4 is computed (tdm26.m:473) but was referenced only in the m<3 branch.
% The m>=4 branch was copied from m=3 and inherited the omission.
%
% WHY IT MATTERS: with only -gam*k4, the gam-dependent force is a pure NEGATIVE
% STIFFNESS, in phase with displacement. A stiffness does zero net work over a
% cycle, so the OHC force could not inject energy no matter how large gam got --
% raising gam could only drive STATIC divergence, preferentially at the
% compliant apex. That is precisely the observed m=3 failure at gam=1 (apical,
% zero latency, +70..+91 dB) and the m=4 collapse at gam=0.5.
%
% pa.r4act = 1 (default, corrected) / 0 (reproduces the old broken behavior).

gv = [0 0.30 0.50 0.70 0.85 1.00];
fp = [1 2 4 8];
ISV = 1391:-10:11;
% INTERIOR ONLY. The dense grid reaches x/L=0.993 (helicotrema) and 0.008
% (stapes), well outside the original saved span 0.057..0.811. In ohc_tip2 the
% peak search landed on those boundary regions and produced lat~-0.01 with
% xbest~0.95 -- a MEASUREMENT artifact that inflated the degeneracy count. Keep
% the dense sampling but score only the interior the original isv covered.
XLO = 0.05; XHI = 0.85;

for m = [3 4]
for ab = [0 1]
    lbl = {'broken control','CORRECTED'};
    fprintf('\n=========== m=%d   r4act=%d  (%s) ===========\n', m, ab, lbl{ab+1});
    fprintf('  gam   contrast(dB)  peak@2k(dB)  gain(dB)  xbest@2k  lat@2k(ms)  degen\n');
    pk0 = NaN;
    for gi = 1:numel(gv)
        pa = modpar26(m); pa.gam = gv(gi); pa.isv = ISV; pa.r4act = ab;
        try
            evalc('S = tdm26(0,pa,0,0);');
        catch e
            fprintf('  %4.2f   FAILED: %s\n', gv(gi), e.message); continue
        end
        nf=numel(S.f); f=S.f(:); P=fft(S.ped); P=P(1:nf);
        np=size(S.d1,2); xp=pa.isv(:)'/pa.n;
        inr = xp>=XLO & xp<=XHI;           % interior mask (see XLO/XHI above)
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
        w=2*pi*f*1000; XB=nan(1,numel(fp)); LT=XB; PKf=XB;
        for k=1:numel(fp)
            [~,ifq]=min(abs(f-fp(k)));
            ex=Am(ifq,:); ex(~isfinite(ex))=0; ex(~inr)=0;   % interior only
            [~,ib]=max(ex); if (ex(ib)<=0), continue; end
            XB(k)=xp(ib); PKf(k)=20*log10(ex(ib));
            ph=unwrap(angle(H(:,ib))); gd=-gradient(ph,w);
            LT(k)=median(gd(max(1,ifq-3):min(nf,ifq+3)))*1000;
        end
        deg = (std(XB,0,'omitnan')<1e-3) || (median(abs(LT),'omitnan')<0.30);
        if (gi==1), pk0=PKf(2); end                  % gam=0 passive reference
        fprintf('  %4.2f   %8.1f     %8.1f   %+7.1f   %7.4f   %8.2f     %d\n', ...
                gv(gi), chi, PKf(2), PKf(2)-pk0, XB(2), LT(2), deg);
    end
end
end
fprintf(['\nPASS for the fix = with r4act=1 the tip GROWS with gam (gain rising,\n' ...
         'contrast rising), latency LENGTHENS and best place moves APICALLY as gam\n' ...
         'rises, and degen stays 0. The broken control should show the opposite or\n' ...
         'a degenerate blow-up.\n']);
disp('R4_TEST_DONE');
