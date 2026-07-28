% Is the compression cliff a PLACE SHIFT?  Track the dominant neural-response
% place vs level.  A smooth compressive nonlinearity moves the place gradually;
% an abrupt basalward jump between two levels IS the cliff.
L=load('parfit26_recip.mat'); pa=L.R.pa; pa.hbmode='bm';
for fr=[0.5 2]
    fprintf('\n%.1f kHz  (3-chamber)\n',fr);
    fprintf('  lvl  lat(ms)   dlat%%   pkplace  pkfrac  centroid  d2/hbmx    ratio\n');
    prev=NaN; pw=NaN;
    for lv=[20 40 60 80]
        p.fr=fr; p.lv=lv; p.pa=pa;
        S=tdm26('wnr1',p,0,0); d=S.dgn; w=max(S.wnr);
        dl=NaN; if (isfinite(prev)), dl=100*(S.tpk-prev)/prev; end
        rt=NaN; if (isfinite(pw)&&pw>0), rt=w/pw; end
        fprintf('  %3.0f %7.2f %7.1f   %6d  %6.3f  %8.3f %8.2f %8.2f\n', ...
                lv,S.tpk,dl,d.pkplace,d.pkfrac,d.pkcen,d.ratio,rt);
        prev=S.tpk; pw=w;
    end
end
fprintf('\npkfrac jumping basalward (toward 0) between two levels => cliff is a PLACE SHIFT\n');
disp('CLIFF_PLACE_DONE');
