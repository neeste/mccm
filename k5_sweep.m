% sgn=+1 INJECTS energy (ohcP>0) but far too weakly to amplify. The OC-height
% impedance k5 sets the PHASE of the OHC force reaching the mechanics, so sweep
% it and find where injected power PEAKS -- that locates the correct phase
% instead of leaving k5 as an unexamined placeholder in the amplifier path.
base=modpar26(4); k5ref=base.k5o;
fprintf('\n  k5/ref |  lat(ms)  d2/hbmx     max|WNR|         ohcP     ratio-vs-passive\n');
fprintf('%s\n',repmat('-',1,78));
% passive reference at same condition
q.fr=2; q.lv=60; q.pa=modpar26(4); q.pa.hbmode='bm'; q.pa.ohcgain=0;
P=tdm26('wnr1',q,0,0); wpass=max(P.wnr);
fprintf('  passive |  %6.2f %8.2f  %11.3e  %+11.3e        1.00\n', P.tpk, P.dgn.ratio, wpass, P.dgn.ohcP);
for sc=[0.01 0.1 0.3 1 3 10 100]
    p.fr=2; p.lv=60; p.pa=modpar26(4); p.pa.hbmode='bm';
    p.pa.ohcsgn=+1; p.pa.ohcgain=1.0; p.pa.k5o=k5ref*sc;
    try
        S=tdm26('wnr1',p,0,0); d=S.dgn; w=max(S.wnr);
        fprintf('  %6.2f |  %6.2f %8.2f  %11.3e  %+11.3e        %.2f\n', ...
                sc, S.tpk, d.ratio, w, d.ohcP, w/wpass);
    catch e
        fprintf('  %6.2f |  FAILED: %s\n', sc, e.message);
    end
end
fprintf('\npeak ohcP => correct OC-height phase.  ratio-vs-passive >1 => amplifying.\n');
disp('K5_SWEEP_DONE');
