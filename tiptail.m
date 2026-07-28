% Tip/tail structure and the delay it produces, from the SAME linearized impulse
% responses.  Q10 conflates tip and tail whenever tip-tail contrast < 10 dB, so
% measure them separately -- and compute the group delay, since the forward
% latency accumulates in the TIP.
%   contrast = |H(BF)| / |H(BF/4)|      (tip above tail, dB)
%   Q3       = BF / (3 dB bandwidth)    (tip sharpness proper)
%   tauBF    = group delay at BF (ms)   -> N = tauBF*BF = delay in CYCLES
% d<1 requires N to RISE with CF.
for nch=[2 3 4]
    fprintf('\n===== nch=%d =====\n',nch);
    try
        pa=modpar26(nch); evalc('S=tdm26(0,nch,0,0);');
        nf=numel(S.f); f=S.f(:); P=fft(S.ped); P=P(1:nf);
        fprintf('  place    BF     contrast   Q3     Q10    tau(ms)     N=tau*BF\n');
        for i=1:size(S.d1,2)
            D=fft(S.d1(:,i)-S.d2(:,i)); D=D(1:nf);
            H=D./max(abs(P),eps); A=abs(H).*f(:);
            ok=f>0.15 & f<18; A(~ok)=0;
            [pk,ip]=max(A); if (pk<=0), continue; end
            bf=f(ip);
            if (bf<0.6 || bf>16), continue; end          % skip range-limit artifacts
            % tail 2 octaves below BF
            [~,it]=min(abs(f-bf/4)); tail=A(it);
            con=20*log10(pk/max(tail,eps));
            bw=@(fr) deal_bw(A,f,ip,pk/10^(fr/20));
            q3=bf/max(bw(3),eps); q10=bf/max(bw(10),eps);
            ph=unwrap(angle(H)); w=2*pi*f*1000;
            gd=-gradient(ph,w);                           % s
            k1=max(1,ip-3); k2=min(nf,ip+3);
            tau=median(gd(k1:k2))*1000;                   % ms
            fprintf('  %5d %7.3f  %8.1f %6.2f %6.2f %9.3f %10.2f\n', ...
                    pa.isv(i), bf, con, q3, q10, tau, tau*bf);
        end
    catch e
        fprintf('  FAILED: %s\n', e.message);
    end
end
fprintf('\ncontrast>>10 dB => distinct tip+tail (Q10 then valid).\n');
fprintf('N rising with BF => the scaling violation that gives d<1.\n');
disp('TIPTAIL_DONE');
