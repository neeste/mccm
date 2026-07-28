% The cliff is a TUNING-BROADENING failure: the model's tuning holds, then gives
% way. Damping sets how near criticality the loop sits, and near criticality the
% response is hypersensitive to gam -> abrupt collapse. More damping should make
% the amplified->unamplified transition GRADUAL, at the cost of amplification.
% Look for a window: smooth steps AND enough total latency change (~-28%/20dB).
% maperr scored too -- smoothing the cliff by wrecking tuning is not a fix.
L=load('parfit26_recip.mat'); pa0=L.R.pa; pa0.hbmode='bm';
fprintf('\nr1 scale |  maperr  | lvl  lat(ms)   dlat%%   centroid  d2/hbmx   ratio\n');
fprintf('%s\n',repmat('-',1,78));
for sc=[0.5 1 2 4]
    pa=pa0; pa.r1o=pa0.r1o*sc;
    try
        Rf=fdm26(struct('pa',pa)); mp=Rf.maperr;
    catch, mp=NaN; end
    prev=NaN; pw=NaN;
    for lv=[20 40 60 80]
        p.fr=2; p.lv=lv; p.pa=pa;
        try
            S=tdm26('wnr1',p,0,0); d=S.dgn; w=max(S.wnr);
            dl=NaN; if (isfinite(prev)), dl=100*(S.tpk-prev)/prev; end
            rt=NaN; if (isfinite(pw)&&pw>0), rt=w/pw; end
            if (lv==20)
                fprintf('  %6.2f | %7.1f  | %3.0f %7.2f %7.1f %10.3f %8.2f %7.2f\n', ...
                        sc,mp,lv,S.tpk,dl,d.pkcen,d.ratio,rt);
            else
                fprintf('         |          | %3.0f %7.2f %7.1f %10.3f %8.2f %7.2f\n', ...
                        lv,S.tpk,dl,d.pkcen,d.ratio,rt);
            end
            prev=S.tpk; pw=w;
        catch e
            fprintf('         |          | %3.0f  FAILED: %s\n',lv,e.message);
        end
    end
end
fprintf('\nSuccess = dlat%% near -28 at EVERY step with maperr still acceptable (<=185).\n');
disp('R1_SWEEP_DONE');
