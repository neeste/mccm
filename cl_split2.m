% PHASE 2, retried on the FIXED topology (d1: ST<->SV, CL a side compartment).
% ohcsplit = fraction of OHC somatic force routed to the cortilymph pump (d3);
% the remainder drives the BM. ohcsplit=0 should now behave much like the
% 3-chamber (genuine regression test); raising it diverts amplifier gain into
% cortilymph pumping.  3-ch ref @2kHz: lat 6.44/4.72/4.40/2.16, ratios 3.08/2.91/10.86
lv=[20 40 60 80];
for fsp=[0 0.5 1.0]
    fprintf('\nohcsplit=%.2f  @2kHz\n',fsp);
    fprintf('  lvl  lat(ms)  d2/hbmx     max|WNR|  ratio\n');
    prev=NaN;
    for j=1:numel(lv)
        p.fr=2; p.lv=lv(j); p.pa=modpar26(4); p.pa.hbmode='bm'; p.pa.ohcsplit=fsp;
        try
            S=tdm26('wnr1',p,0,0); d=S.dgn; mx=max(S.wnr);
            r=NaN; if (isfinite(prev)&&prev>0), r=mx/prev; end
            fprintf('  %3.0f %7.2f %8.2f  %11.3e %6.2f\n', lv(j), S.tpk, d.ratio, mx, r);
            prev=mx;
        catch e
            fprintf('  %3.0f  FAILED: %s\n', lv(j), e.message);
        end
    end
end
fprintf('\nratio<10 => compression active. Does lower ohcsplit restore 3-chamber-like gain?\n');
disp('CL_SPLIT2_DONE');
