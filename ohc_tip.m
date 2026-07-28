% OHC FORCE AMPLITUDE -> TIP, BEST PLACE, AND LATENCY  (m = 2, 3, 4)
%
% Premise (SN): the frequency-place map is LEVEL DEPENDENT -- the maximum of the
% excitation pattern shifts with level.  The underlying cause is the OHC force:
% as its amplitude falls, the tip should collapse and BOTH the latency and the
% best place should shift.  Classic direction: peak moves BASALWARD (isv runs
% apical->basal, so toward HIGHER index) and latency SHORTENS.
%
% KNOB: pa.gam.  cp.gm = pa.gam*pa.gampro (line 501) feeds gam, which enters as
% -gam.*cp.k4 in k_act in BOTH the m<4 branch (654) and the m>=4 branch (677).
% It is therefore the SAME amplifying element in all three models.  pa.ohcgain
% is NOT usable here: it exists only for m>=4 and scales the whole force pair
% including the passive gh.*k3 coupling, confounding gain with coupling.
%
% modpar26 leaves hbnl=0, so the click path is LINEAR: gam is not further
% reduced by the hbt compression of lines 644-645.  Amplifier gain is thus
% varied independently of level, which is the point -- we want the CAUSE of the
% level dependence, not the level dependence itself.
%
% Measured on d1 (BM), consistent with the BM-peak detector fix.

gv  = [1.00 0.70 0.50 0.30 0.15 0.00];    % OHC force amplitude
fp  = [1 2 4 8];                          % probe frequencies (kHz)
mv  = [2 3 4];

R = struct();
for mi = 1:numel(mv)
    m = mv(mi);
    fprintf('\n================ m = %d ================\n', m);
    XB=nan(numel(gv),numel(fp)); LT=XB; PK=XB; BFxb=XB; CHI=nan(numel(gv),1); DD=CHI;
    bfmap=[]; xpos=[];
    for gi = 1:numel(gv)
        pa = modpar26(m); pa.gam = gv(gi);
        try
            evalc('S = tdm26(0,pa,0,0);');
        catch e
            fprintf('  gam=%.2f  FAILED: %s\n', gv(gi), e.message); continue
        end
        nf=numel(S.f); f=S.f(:); P=fft(S.ped); P=P(1:nf);
        np=size(S.d1,2);
        if (isempty(xpos)), xpos = pa.isv(:)'/pa.n; end          % normalized place
        H=zeros(nf,np);
        for i=1:np
            D=fft(S.d1(:,i)); H(:,i)=D(1:nf)./max(abs(P),eps);   % BM transfer
        end
        A=abs(H);
        % ---- per-place BF map (only needed once, at full gain, as the ruler) --
        bf=nan(1,np); co=nan(1,np);
        for i=1:np
            a=A(:,i).*f; a(~(f>0.15&f<18))=0;
            [pk,ip]=max(a); if (pk<=0), continue; end
            if (f(ip)<0.6||f(ip)>16), continue; end
            bf(i)=f(ip);
            [~,it]=min(abs(f-bf(i)/4)); co(i)=20*log10(pk/max(a(it),eps));
        end
        if (gi==1), bfmap=bf; end
        ok=isfinite(co); if (any(ok)), CHI(gi)=co(find(ok,1,'last')); end
        % ---- excitation pattern: best place and latency at each probe freq ----
        w=2*pi*f*1000;
        for k=1:numel(fp)
            [~,ifq]=min(abs(f-fp(k)));
            ex=A(ifq,:);                       % excitation pattern along place
            ex(~isfinite(ex))=0;
            [~,ib]=max(ex); if (ex(ib)<=0), continue; end
            XB(gi,k)=xpos(ib); PK(gi,k)=20*log10(ex(ib));
            if (~isempty(bfmap) && isfinite(bfmap(ib))), BFxb(gi,k)=bfmap(ib); end
            ph=unwrap(angle(H(:,ib))); gd=-gradient(ph,w);
            LT(gi,k)=median(gd(max(1,ifq-3):min(nf,ifq+3)))*1000;   % ms
        end
        fprintf('  gam=%.2f done\n', gv(gi));
    end
    R(mi).m=m; R(mi).XB=XB; R(mi).LT=LT; R(mi).PK=PK; R(mi).BFxb=BFxb; R(mi).CHI=CHI;

    % ---------------- report ----------------
    for k=1:numel(fp)
        fprintf('\n  --- m=%d, probe %.0f kHz ---\n', m, fp(k));
        fprintf('   gam   xbest(x/L)  shift(oct)   lat(ms)   dlat(ms)   peak(dB)  gain(dB)\n');
        for gi=1:numel(gv)
            so=NaN;
            if (isfinite(BFxb(gi,k))&&isfinite(BFxb(1,k))), so=log2(BFxb(gi,k)/BFxb(1,k)); end
            dl=LT(gi,k)-LT(1,k);
            gn=PK(gi,k)-PK(end,k);      % re passive (gam=0)
            fprintf('  %4.2f    %6.3f     %+6.2f     %6.2f    %+6.2f     %6.1f    %+6.1f\n', ...
                    gv(gi), XB(gi,k), so, LT(gi,k), dl, PK(gi,k), gn);
        end
    end
    fprintf('\n  tip-tail contrast at highest CF (dB) vs gam:');
    fprintf(' %.1f', CHI); fprintf('\n');
end
save('ohc_tip.mat','R','gv','fp','mv');
disp('OHC_TIP_DONE');
