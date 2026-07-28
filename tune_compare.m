% Tuning curves from LINEARIZED impulse (short-chirp) responses, 2/3/4 chambers.
% H_i(f) = FFT(d1-d2)_i / FFT(ped)  -- hair-bundle response re eardrum pressure.
% Computed directly rather than via the built-in bf/Qe, which fails at the
% array ends. Q10 = BF / (bandwidth 10 dB below peak).
res=struct();
for nch=[2 3 4]
    fprintf('\n===== nch=%d =====\n',nch);
    try
        pa=modpar26(nch);
        evalc('S=tdm26(0,nch,0,0);');
        nt=size(S.d1,1); nf=numel(S.f); f=S.f(:);
        P=fft(S.ped); P=P(1:nf);
        fprintf('  place   BF(kHz)   peak(dB)    Q10     (parlab %s)\n', pa.parlab);
        BF=nan(1,7); PK=nan(1,7); Q=nan(1,7);
        for i=1:size(S.d1,2)
            D=fft(S.d1(:,i)-S.d2(:,i)); D=D(1:nf);
            H=abs(D./max(abs(P),eps));
            v=H(:).*f(:);                       % velocity-like
            ok=f>0.1 & f<20;                    % ignore DC and the Nyquist end
            vv=v; vv(~ok)=0;
            [pk,ip]=max(vv);
            if (pk<=0), fprintf('  %5d      (no peak)\n', pa.isv(i)); continue; end
            bf=f(ip); thr=pk/sqrt(10);          % 10 dB down
            lo=ip; while (lo>1 && vv(lo)>thr), lo=lo-1; end
            hi=ip; while (hi<nf && vv(hi)>thr), hi=hi+1; end
            bw=f(hi)-f(lo); q=NaN; if (bw>0), q=bf/bw; end
            BF(i)=bf; PK(i)=20*log10(pk); Q(i)=q;
            fprintf('  %5d   %7.3f   %8.1f  %6.2f\n', pa.isv(i), bf, 20*log10(pk), q);
        end
        res(nch).BF=BF; res(nch).PK=PK; res(nch).Q=Q; res(nch).isv=pa.isv;
    catch e
        fprintf('  FAILED: %s\n', e.message);
    end
end
fprintf('\n--- summary: median Q10 over places with a valid peak ---\n');
for nch=[2 3 4]
    if (numel(res)>=nch && ~isempty(res(nch).Q))
        q=res(nch).Q; q=q(isfinite(q));
        fprintf('  nch=%d : median Q10 = %.2f   (n=%d places)\n', nch, median(q), numel(q));
    end
end
save('tune_compare.mat','res');
disp('TUNE_COMPARE_DONE');
