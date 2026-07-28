% DOES THE 3-CHAMBER'S FAILURE COME FROM ITS MICROMECHANICS OR ITS CHAMBERS?
%
% m=2 (works: +39 dB tip) drives the active force from the RELATIVE displacement
%   d3 = d1 - d2      (tdm26.m m<3 branch:  s4tmp = cp.k4.*d3 + cp.r4.*v3)
% m=3 (fails) drives it from d2 ALONE
%   (k_act = gh*k3 - gam*k4, entering s1 as +gam*k4*d2)
% Energy injection depends on the PHASE of the active force against BM velocity,
% so the driving coordinate is exactly what decides amplify-vs-dissipate.
%
% pa.m3form=1 runs the m=2-style micromechanics with the 3-chamber hydrodynamics
% UNCHANGED. That isolates micromechanics from chamber count -- which the earlier
% swap test could not do, since it varied chamber count and parameter set at once.
%
% m=2 is included as the reference the m3form=1 arm should move toward.
% Scored on the interior only (x/L in [0.05,0.85]): the dense grid reaches the
% helicotrema and stapes, where the peak search picks up boundary artifacts.

gv  = [0 0.30 0.50 0.70 0.85 1.00];
ISV = 1391:-10:11; XLO=0.05; XHI=0.85;
cfg = { 'm=2 reference     ', 2, 0
        'm=3 d2-only (curr)', 3, 0
        'm=3 d3=d1-d2 (new)', 3, 1 };

for c = 1:size(cfg,1)
    nm = cfg{c,1}; m = cfg{c,2}; mf = cfg{c,3};
    fprintf('\n===== %s  (m=%d, m3form=%d) =====\n', nm, m, mf);
    fprintf('  gam   contrast(dB)  peak@2k(dB)  gain(dB)  xbest@2k  lat@2k(ms)  degen\n');
    pk0 = NaN;
    for gi = 1:numel(gv)
        pa = modpar26(m); pa.gam = gv(gi); pa.isv = ISV;
        if (m==3), pa.m3form = mf; end
        try
            evalc('S = tdm26(0,pa,0,0);');
        catch e
            fprintf('  %4.2f   FAILED: %s\n', gv(gi), e.message); continue
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
        w=2*pi*f*1000; XB=nan(1,4); LT=XB; PKf=XB; fp=[1 2 4 8];
        for k=1:4
            [~,ifq]=min(abs(f-fp(k)));
            ex=Am(ifq,:); ex(~isfinite(ex))=0; ex(~inr)=0;
            [~,ib]=max(ex); if (ex(ib)<=0), continue; end
            XB(k)=xp(ib); PKf(k)=20*log10(ex(ib));
            ph=unwrap(angle(H(:,ib))); gd=-gradient(ph,w);
            LT(k)=median(gd(max(1,ifq-3):min(nf,ifq+3)))*1000;
        end
        deg = (std(XB,0,'omitnan')<1e-3) || (median(abs(LT),'omitnan')<0.30);
        if (gi==1), pk0=PKf(2); end
        fprintf('  %4.2f   %8.1f     %8.1f   %+7.1f   %7.4f   %8.2f     %d\n', ...
                gv(gi), chi, PKf(2), PKf(2)-pk0, XB(2), LT(2), deg);
    end
end
fprintf(['\nHYPOTHESIS CONFIRMED if the m3form=1 arm gains a real tip: gain RISING\n' ...
         'with gam (m=2 reaches ~+39 dB), contrast RISING, latency LENGTHENING,\n' ...
         'best place moving APICALLY, degen 0 throughout -- i.e. it moves toward\n' ...
         'the m=2 reference rather than tracking the d2-only arm.\n' ...
         'REFUTED if it still shows falling contrast / the apical zero-latency\n' ...
         'blow-up, which would put the failure in the CHAMBERS, not the force law.\n']);
disp('M3FORM_TEST_DONE');
