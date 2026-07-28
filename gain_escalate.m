% k5 tuning gave ~1.2x amplification at any phase. The remaining untested axis:
% the INJECTING sign (sgn=+1) has only ever been run at gain<=1.0. Injected power
% scales with ohcgain, so push it. (The gain=1.0 runaway was sgn=-1, the
% DISSIPATING sense, so it does not bound the injecting sign.)
% Look for: ohcP rising, WNR/passive rising well above 1.3, ratio falling below 10.
k5s=0.3;   % near the ohcP peak
q.fr=2; q.lv=60; q.pa=modpar26(4); q.pa.hbmode='bm'; q.pa.ohcgain=0;
P=tdm26('wnr1',q,0,0); wpass=max(P.wnr);
fprintf('\n passive: WNR=%.3e  d2/hbmx=%.2f  lat=%.2f\n', wpass, P.dgn.ratio, P.tpk);
fprintf('\n  gain |  lat(ms)   d2/hbmx      max|WNR|         ohcP    vs-passive\n');
fprintf('%s\n',repmat('-',1,72));
for g=[1 3 10 30 100]
    p.fr=2; p.lv=60; p.pa=modpar26(4); p.pa.hbmode='bm';
    p.pa.ohcsgn=+1; p.pa.ohcgain=g; p.pa.k5o=modpar26(4).k5o*k5s;
    try
        S=tdm26('wnr1',p,0,0); d=S.dgn; w=max(S.wnr);
        fprintf('  %5.0f | %7.2f %10.2f  %11.3e  %+11.3e   %9.2f\n', ...
                g, S.tpk, d.ratio, w, d.ohcP, w/wpass);
    catch e
        fprintf('  %5.0f | FAILED: %s\n', g, e.message);
    end
end
fprintf('\nIf WNR/passive stays ~1.2 while ohcP scales, the pump injects energy but\n');
fprintf('cannot convert it into BM amplification => negative result for the mechanism.\n');
disp('GAIN_ESCALATE_DONE');
