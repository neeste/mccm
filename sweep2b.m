% Rerun nch=2 with maptol at its own frontier (105, not 167), matching the nch=1
% rerun that landed maperr 90.97 -- a new project best, beating 104.63.
% nch=2 is the same model as nch=1 by the symmetry constraint, so this both
% corrects the fit AND extends the m=1 == m=2 identity result to the corrected
% one. The prior pair landed identically at 129.96; if these land identically at
% ~90.97 the identity survives a properly-constrained fit.
t0=tic;
try
    R = parfit26(2, struct('maxfe',150,'wgain',0.01,'out','sweep_nch2b.mat'));
    S = score26(R.pa,'fast',false); m = abr_metric(R.pa,false);
    fprintf('\n  DONE nch=2 | maperr %.2f | amp_gain %+.2f | slope %.3f | maxRe_osc %.1f | %.1f h\n', ...
            S.maperr, S.amp_gain, m.slope, S.maxRe_osc, toc(t0)/3600);
    fprintf('    J terms: slope %.4f  map %.4f  gain %.4f  osc %.4f\n', ...
            abs(m.slope-0.413), 0.001*max(0,S.maperr-105), ...
            0.01*max(0,40-S.amp_gain), 0.005*max(0,S.maxRe_osc+40));
    fprintf('    nch=1 rerun for comparison: maperr 90.97, amp +54.54, slope 0.606\n');
    fprintf('    prior nch=2 (maptol 167):   maperr 129.96, slope 0.596\n');
catch e
    fprintf('  FAILED after %.1f h: %s\n', toc(t0)/3600, e.message(1:min(100,end)));
end
