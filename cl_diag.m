% Does the 4-chamber actually have a TRAVELING WAVE, or is the short 1.90 ms
% latency a basal artifact?  If a wave propagates to a proper CF place, latency
% must FALL systematically with frequency (as it does in the 3-chamber).
% Also check stability, and compare WNR magnitudes 3-ch vs 4-ch.
fprintf('\n f(kHz) |  3-chamber lat / max      |  4-chamber lat / max\n');
fprintf('%s\n',repmat('-',1,62));
L=load('parfit26_recip.mat'); pa3=L.R.pa; pa3.hbmode='bm';
pa4=modpar26(4); pa4.hbmode='bm';
for fr=[0.5 1 2 4]
    a.fr=fr; a.lv=60; a.pa=pa3; S3=tdm26('wnr1',a,0,0);
    b.fr=fr; b.lv=60; b.pa=pa4;
    try
        S4=tdm26('wnr1',b,0,0);
        fprintf('%6.2f  |  %6.2f ms  %9.2e    |  %6.2f ms  %9.2e\n', ...
                fr, S3.tpk, max(S3.wnr), S4.tpk, max(S4.wnr));
    catch e
        fprintf('%6.2f  |  %6.2f ms  %9.2e    |  FAILED: %s\n', fr, S3.tpk, max(S3.wnr), e.message);
    end
end
fprintf('\n(3-chamber latency falls 0.5->4 kHz; the 4-chamber must too if a wave propagates)\n');
fprintf('\n--- 4-chamber stability ---\n');
try
    evalc('S=tdm26(''coupeig'',struct(''pa'',pa4));');
    fprintf('coupeig maxRe_osc = %+.1f  (%s)\n', S.maxRe_osc, ...
            char(74*(S.maxRe_osc<0)+85*(S.maxRe_osc>=0)));
    if (S.maxRe_osc<0), fprintf('  sub-critical (stable)\n'); else, fprintf('  UNSTABLE\n'); end
catch e
    fprintf('coupeig failed: %s\n', e.message);
end
disp('CL_DIAG_DONE');
