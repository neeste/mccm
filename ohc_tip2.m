% OHC FORCE -> MAP SHIFT, refined.   Follow-up to ohc_tip.m.
%
% Two limitations of the first pass, both fixed here:
%  (1) pa.isv saved only 7 places (x/L 0.811..0.057, ~1 octave apart), so every
%      measured shift was exactly one index step -- the minimum resolvable.
%      isv is only a RECORDING choice (par_CEL16.m:14 "BM locations to save"),
%      not physics, so a dense grid changes nothing about the model.
%  (2) gam was sampled too coarsely to see inside the cliff: for m=2 almost the
%      entire tip collapse happened between gam=1.0 and gam=0.7.
%
% MAP SHIFT is now measured properly.  Rather than tracking one best place, fit
% the WHOLE place-frequency map at each gain:  log2(BF) = a + b*(x/L).  A real
% level/gain-dependent map shows up as a change in the intercept a, in octaves,
% using all valid places at once.  That is immune to place quantization.
%
% Degeneracy flag replaces coupeig here: coupeig would need eig of a 2*n*dof
% system (8406 for m=4), far too slow for 28 runs.  The degenerate cases in the
% first pass had an unmistakable signature -- best place identical at every
% probe frequency AND |latency| ~ 0 -- so that is what is flagged, and coupeig
% can be run afterwards on the few cases that matter.

cfg = { 2, [1.00 0.98 0.95 0.92 0.90 0.85 0.80 0.70 0.50 0.20 0.00]
        3, [1.00 0.95 0.90 0.80 0.70 0.50 0.30 0.00]
        4, [1.00 0.95 0.90 0.85 0.80 0.75 0.70 0.60 0.50] };
fp = [1 2 4 8];
ISV = 1391:-10:11;                 % 139 places, apical -> basal (as isv runs)

R = struct();
for ci = 1:size(cfg,1)
    m = cfg{ci,1}; gv = cfg{ci,2}; ng = numel(gv);
    fprintf('\n================ m = %d  (%d places, %d gains) ================\n', ...
            m, numel(ISV), ng);
    A0=nan(ng,1); B0=nan(ng,1); R2=nan(ng,1); CHI=nan(ng,1);
    XB=nan(ng,numel(fp)); LT=XB; PK=XB; DEG=false(ng,1);
    for gi = 1:ng
        pa = modpar26(m); pa.gam = gv(gi); pa.isv = ISV;
        try
            evalc('S = tdm26(0,pa,0,0);');
        catch e
            fprintf('  gam=%.2f FAILED: %s\n', gv(gi), e.message); continue
        end
        nf=numel(S.f); f=S.f(:); P=fft(S.ped); P=P(1:nf);
        np=size(S.d1,2); xp=pa.isv(:)'/pa.n;
        H=zeros(nf,np);
        for i=1:np, D=fft(S.d1(:,i)); H(:,i)=D(1:nf)./max(abs(P),eps); end
        Am=abs(H);
        % ---- dense BF map, then fit log2(BF) = a + b*(x/L) -------------------
        bf=nan(1,np); co=nan(1,np);
        for i=1:np
            a=Am(:,i).*f; a(~(f>0.15&f<18))=0;
            [pk,ip]=max(a); if (pk<=0), continue; end
            if (f(ip)<0.5||f(ip)>16.5), continue; end
            bf(i)=f(ip); [~,it]=min(abs(f-bf(i)/4));
            co(i)=20*log10(pk/max(a(it),eps));
        end
        ok=isfinite(bf);
        if (sum(ok)>=10)
            pf=polyfit(xp(ok),log2(bf(ok)),1); B0(gi)=pf(1); A0(gi)=pf(2);
            rr=log2(bf(ok))-polyval(pf,xp(ok));
            R2(gi)=1-sum(rr.^2)/max(sum((log2(bf(ok))-mean(log2(bf(ok)))).^2),eps);
            iv=find(ok); CHI(gi)=co(iv(end));       % contrast at most basal valid
        end
        % ---- best place (parabolic-interpolated) and latency per probe freq ---
        w=2*pi*f*1000;
        for k=1:numel(fp)
            [~,ifq]=min(abs(f-fp(k)));
            ex=Am(ifq,:); ex(~isfinite(ex))=0;
            [~,ib]=max(ex); if (ex(ib)<=0), continue; end
            xb=xp(ib);
            if (ib>1 && ib<np)                       % parabolic peak refinement
                y1=log(max(ex(ib-1),eps)); y2=log(max(ex(ib),eps)); y3=log(max(ex(ib+1),eps));
                den=(y1-2*y2+y3); if (abs(den)>eps)
                    dl=0.5*(y1-y3)/den; dl=max(min(dl,1),-1);
                    xb=xp(ib)+dl*(xp(min(ib+1,np))-xp(max(ib-1,1)))/2;
                end
            end
            XB(gi,k)=xb; PK(gi,k)=20*log10(ex(ib));
            ph=unwrap(angle(H(:,ib))); gd=-gradient(ph,w);
            LT(gi,k)=median(gd(max(1,ifq-3):min(nf,ifq+3)))*1000;
        end
        DEG(gi) = (std(XB(gi,:),0,'omitnan')<1e-3) || (median(abs(LT(gi,:)),'omitnan')<0.30);
        fprintf('  gam=%.2f done%s\n', gv(gi), char(32*(1-DEG(gi))*ones(1,~DEG(gi))) );
        if (DEG(gi)), fprintf('      ^ DEGENERATE (place freq-independent or lat~0)\n'); end
    end
    R(ci).m=m; R(ci).gv=gv; R(ci).A0=A0; R(ci).B0=B0; R(ci).R2=R2;
    R(ci).CHI=CHI; R(ci).XB=XB; R(ci).LT=LT; R(ci).PK=PK; R(ci).DEG=DEG;

    fprintf('\n  --- m=%d: place-frequency MAP vs OHC force ---\n', m);
    fprintf('   gam   mapshift(oct)  slope b   R2    contrast(dB)  degen\n');
    for gi=1:ng
        fprintf('  %4.2f     %+7.3f     %+6.2f  %5.2f    %7.1f      %d\n', ...
                gv(gi), A0(gi)-A0(1), B0(gi), R2(gi), CHI(gi), DEG(gi));
    end
    for k=1:numel(fp)
        fprintf('\n  --- m=%d, probe %.0f kHz ---\n', m, fp(k));
        fprintf('   gam   xbest(x/L)  dx(basal+)   lat(ms)   dlat(ms)   peak(dB)\n');
        for gi=1:ng
            fprintf('  %4.2f    %6.4f     %+7.4f    %6.2f    %+6.2f    %7.1f\n', ...
                    gv(gi), XB(gi,k), XB(1,k)-XB(gi,k), LT(gi,k), LT(gi,k)-LT(1,k), PK(gi,k));
        end
    end
end
save('ohc_tip2.mat','R','fp','ISV');
disp('OHC_TIP2_DONE');
