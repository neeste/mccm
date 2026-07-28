% hbsc sets WHERE the amplifier disengages, not how abruptly. Abruptness comes
% from a near-critical loop -- and every place disengages TOGETHER because
% pa.gampro is uniform ones(n,1) (never fitted). A GRADED gain profile should
% stagger saturation across place and smooth the aggregate transition.
%   gampro = exp(g*(xf-0.5)):  g>0 apical-biased, g<0 basal-biased, g=0 uniform
L=load('parfit26_recip.mat'); pa0=L.R.pa; pa0.hbmode='bm';
n=pa0.n; xf=((0:n-1)')/(n-1);
for g=[0 1 2 -1 -2]
    fprintf('\ngampro grade g=%+.0f   (2 kHz)   [target: uniform -28%%/20dB]\n',g);
    fprintf('  lvl  lat(ms)   dlat%%   centroid  d2/hbmx   ratio\n');
    prev=NaN; pw=NaN;
    for lv=[20 40 60 80]
        p.fr=2; p.lv=lv; p.pa=pa0; p.pa.gampro=exp(g*(xf-0.5));
        try
            S=tdm26('wnr1',p,0,0); d=S.dgn; w=max(S.wnr);
            dl=NaN; if (isfinite(prev)), dl=100*(S.tpk-prev)/prev; end
            rt=NaN; if (isfinite(pw)&&pw>0), rt=w/pw; end
            fprintf('  %3.0f %7.2f %7.1f %10.3f %8.2f %7.2f\n',lv,S.tpk,dl,d.pkcen,d.ratio,rt);
            prev=S.tpk; pw=w;
        catch e
            fprintf('  %3.0f  FAILED: %s\n',lv,e.message);
        end
    end
end
fprintf('\nSuccess = dlat%% closer to -28 at EVERY step (staggered disengagement).\n');
fprintf('CAVEAT: grading gampro also alters tuning, so maperr must be rechecked.\n');
disp('GAMPRO_SWEEP_DONE');
