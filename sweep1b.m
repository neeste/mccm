% Rerun nch=1 with maptol at its own frontier (105, not 167). The previous run
% had map=0.0000 for all 150 evals and let maperr drift 104.6 -> 129.96 free.
t0=tic;
try
    R = parfit26(1, struct('maxfe',150,'wgain',0.01,'out','sweep_nch1b.mat'));
    S = score26(R.pa,'fast',false); m = abr_metric(R.pa,false);
    fprintf('\n  DONE nch=1 | maperr %.2f | amp_gain %+.2f | slope %.3f | maxRe_osc %.1f | %.1f h\n', ...
            S.maperr, S.amp_gain, m.slope, S.maxRe_osc, toc(t0)/3600);
    fprintf('    J terms: slope %.4f  map %.4f  gain %.4f  osc %.4f\n', ...
            abs(m.slope-0.413), 0.001*max(0,S.maperr-105), ...
            0.01*max(0,40-S.amp_gain), 0.005*max(0,S.maxRe_osc+40));
    fprintf('    prior run (maptol 167): maperr 129.96, slope 0.596, map term 0.0000\n');
catch e
    fprintf('  FAILED after %.1f h: %s\n', toc(t0)/3600, e.message(1:min(100,end)));
end
