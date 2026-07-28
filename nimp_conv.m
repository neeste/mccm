% IS nimp THE REASON THE MODEL DOES NOT CONVERGE IN n?
%
% SN's lead. pa.nimp is the corrector-iteration count in tdm26's predictor-
% corrector step (tdm26.m:154-166):
%     [cur,ac] = accel(...,cur,cur);      predictor acceleration
%     nxt = update(pa,cur,nxt,ac,ac);     explicit predict
%     for imp=1:pa.nimp
%         [cur,an] = accel(...,cur,nxt);  re-evaluate at the predicted state
%         nxt = update(pa,cur,nxt,ac,an); trapezoidal correct
%     end
% At nimp=1 (the default everywhere in this project) the implicit trapezoidal
% equation gets ONE correction rather than being iterated to convergence. In a
% stiff near-critical system that leaves a residual which DOES NOT SHRINK when
% dx is refined -- exactly the signature conv_slope found.
%
% WHAT conv_slope FOUND, and why another explanation is needed. Refining n alone
% left both models unconverged, m=3 worst: slope 0.117 -> 0.400 -> 0.255
% (non-monotonic), maperr 401.8 -> 499.3 -> 790.4, level_c 31.96 -> 17.08 ->
% 29.92. I attributed that to dt being held FIXED at 2e-6 while dx halved, which
% varies dt/dx fourfold and is a real confound. But an under-converged corrector
% would produce the same non-convergence INDEPENDENTLY of dt/dx, and it is much
% cheaper to test.
%
% PART A asks whether nimp moves anything at all at fixed n. If the metrics are
% flat in nimp, one correction is already converged and this is not the lever.
% PART B is the real question: does a converged corrector make the model
% CONVERGE IN n? Part A can move the numbers without fixing the n-dependence;
% only Part B settles it.
%
% m=3 IS USED THROUGHOUT because abr_metric costs ~102 s there against ~2290 s
% at m=4, and because m=3 showed the WORST grid behaviour, so any fix must show
% up here first. If nimp helps, m=4 gets the same treatment afterwards.
%
% COST higher nimp means more accel calls per step, so nimp=5 is roughly 3-5x
% the march cost. Still minutes at m=3.

NIMP = [1 2 3 5];
NN   = [701 1401 2801];
NIMP_B = 4;                       % Part B corrector count

fprintf('\n  PART A -- does nimp move the metrics at fixed n=1401?\n');
fprintf('  reference (nimp=1, n=1401): slope 0.400  level_c 17.08  maperr 499.3\n\n');
fprintf('  nimp | slope | level_c | maperr\n');
fprintf('%s\n', repmat('-',1,40));
for k = NIMP
    pa = modpar26(3); pa.nimp = k;
    sl = NaN; lc = NaN; mp = NaN;
    try
        evalc('m = abr_metric(pa,false);'); sl = m.slope; lc = m.level_c;
    catch e
        fprintf('  %4d | abr FAILED: %s\n', k, e.message);
    end
    try, S = score26(pa,'fast',false); mp = S.maperr; catch, end
    fprintf('  %4d | %5.3f | %7.2f | %6.1f\n', k, sl, lc, mp);
end

fprintf('\n  PART B -- with nimp=%d, does the model converge in n?\n', NIMP_B);
fprintf('  at nimp=1 this was: 0.117 / 0.400 / 0.255 (non-monotonic, unconverged)\n\n');
fprintf('  n    | slope | level_c | maperr\n');
fprintf('%s\n', repmat('-',1,40));
for nn = NN
    pa = modpar26(3); pa.nimp = NIMP_B;
    if (nn ~= pa.n)
        try, pa = setn(pa, nn); catch e
            fprintf('  %4d | setn FAILED: %s\n', nn, e.message); continue
        end
    end
    sl = NaN; lc = NaN; mp = NaN;
    try
        evalc('m = abr_metric(pa,false);'); sl = m.slope; lc = m.level_c;
    catch e
        fprintf('  %4d | abr FAILED: %s\n', nn, e.message);
    end
    try, S = score26(pa,'fast',false); mp = S.maperr; catch, end
    fprintf('  %4d | %5.3f | %7.2f | %6.1f\n', nn, sl, lc, mp);
end

fprintf(['\n  PART A flat -> one correction is already converged; nimp is not the\n' ...
         '    lever and the dt/dx confound remains the best explanation.\n' ...
         '  PART A moves, PART B still scattered -> the corrector matters but does\n' ...
         '    not fix the grid dependence; dt must be scaled with dx as well.\n' ...
         '  PART A moves, PART B tight -> SN is right, the model was never\n' ...
         '    converged in TIME, and every metric in this project has been read\n' ...
         '    off an under-solved time step. That would make nimp a prerequisite\n' ...
         '    for any fit, and would explain the unphysical slopes (negative d)\n' ...
         '    seen in vent_slope and row4_slope.\n']);
disp('NIMP_CONV_DONE');
