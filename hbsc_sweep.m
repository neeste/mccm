% The cliff = basal spread at high level, caused by compression FAILING there:
% gam = gm/(1+hbsc*log hbt) gives only 17% reduction at hbt=161 with hbsc=0.04.
% Sweep hbsc: stronger high-level gain reduction should hold the response down,
% prevent the basal spread, and smooth the latency steps toward a uniform -28%.
L=load('parfit26_recip.mat'); pa0=L.R.pa; pa0.hbmode='bm';
for hb=[0.04 0.1 0.3 1.0]
    fprintf('\nhbsc=%.2f   (2 kHz)   [target: uniform -28%%/20dB]\n',hb);
    fprintf('  lvl  lat(ms)   dlat%%   centroid  d2/hbmx   ratio   gam-red%%\n');
    prev=NaN; pw=NaN;
    for lv=[20 40 60 80]
        p.fr=2; p.lv=lv; p.pa=pa0; p.pa.hbsc=hb;
        try
            S=tdm26('wnr1',p,0,0); d=S.dgn; w=max(S.wnr);
            dl=NaN; if (isfinite(prev)), dl=100*(S.tpk-prev)/prev; end
            rt=NaN; if (isfinite(pw)&&pw>0), rt=w/pw; end
            gr=100*(1-1/(1+hb*log(max(d.ratio,1))));
            fprintf('  %3.0f %7.2f %7.1f %10.3f %8.2f %7.2f %9.1f\n', ...
                    lv,S.tpk,dl,d.pkcen,d.ratio,rt,gr);
            prev=S.tpk; pw=w;
        catch e
            fprintf('  %3.0f  FAILED: %s\n',lv,e.message);
        end
    end
end
fprintf('\nSuccess = dlat%% converging on -28 at every step, centroid moving smoothly.\n');
disp('HBSC_SWEEP_DONE');
