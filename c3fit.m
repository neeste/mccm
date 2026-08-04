% C3FIT -- the Banff-window run. nch=3, 17 params, hbsc free for the first time.
%
% ####################################################################
% ### CORRECTION 2026-08-03, AFTER THE RUN. THE PREMISE BELOW IS   ###
% ### WRONG. THE LEVER SIGNS DO NOT INVERT.                        ###
% ####################################################################
%
% lcprobe.m swept hbsc at nch=3 from the finished fit (anchor reproduced
% seg2 exactly; maperr/osc/amp verified CONSTANT to 0.00e+00 across the
% sweep, since neither hbsc nor hbmx reaches fdm26 or coupeig):
%
%     hbsc     slope   lvl_c        hbsc     slope   lvl_c
%     0.005   0.6506   1.802        0.06    0.3665  12.499
%     0.010   0.6380   2.403        0.10    0.4379  12.230
%     0.020   0.5903   4.189        0.20    0.4452  14.907
%     0.030   0.5107   5.632
%     0.0406  0.4151   6.469  <- seg2
%
% Over 0.005 -> 0.06 hbsc LOWERS slope and RAISES level_c -- the SAME
% directions as at nch=1. The axis is non-monotonic with a minimum near
% hbsc 0.06 (slope 0.3665). c3probe sampled only 0.04 and 0.65, ONE ON
% EACH SIDE of that minimum, saw slope rise, and read it as a sign flip.
% It is a sampling artifact of a non-monotonic axis, not a structural
% difference between chamber counts.
%
% Same failure class as the +240 maxRe_osc figure: a two-point reading
% promoted to a mechanism claim without a sweep over the range the claim
% needs. abr-tuning-levers already carries the rule -- "ALWAYS verify an
% atlas row with a direct sweep over the range the claim requires."
%
% ALSO WRONG IN PRACTICE: this run was built to free hbsc, and hbsc never
% moved. Across all three segments it went 0.0404 / 0.0406 / 0.0407 from
% a 0.0400 default, and ace moved 0.83%. No parameter moved more than
% 3.71% (r3o). J fell 0.1152 -> 0.0012 on a sub-4% collective nudge of
% the ordinary impedance vector; the stated premise was never exercised.
%
% WHAT THE FIT ACTUALLY SPENT, none of it charged by J: level_c 4.796 ->
% 6.399 (1.580 -> 1.869 %/dB against a 1.62 target) and oscillatory
% margin -220.3 -> -53.7 dB. Both are threshold penalties sitting in
% their dead zones, so J 0.0012 reads "solved" largely because three of
% five terms are inside tolerance bands rather than good. seg2 at
% level_c 6.469 sits just under a CLIFF -- level_c jumps to 12.499 by
% hbsc 0.060.
% ####################################################################
%
% ============================ WHY THIS SHAPE ============================
%
% [SUPERSEDED -- see the correction above. Kept because it is why the run
% was shaped this way, not because it is true.]
%
% THE nch=1 LEVERS DO NOT TRANSFER, measured 2026-08-02 (c3probe). Transplanting
% ace and hbsc from the nch=1 optimum onto the best nch=3 fit makes every version
% WORSE, and the signs INVERT:
%
%   config              slope   lvl_c   maperr    amp     osc   |    J
%   nch=3 as-is        0.5023   4.796   147.05  +37.42  -220.3  | 0.1152
%   + hbsc only        0.5961  13.865   147.05  +37.42  -220.3  | 0.9454
%   + ace only         0.6253   2.411   199.18  +28.82  -253.0  | 0.4822
%   + both             0.6062   1.665   199.18  +28.82  -253.0  | 0.5377
%
% At nch=1, ace LOWERED the slope and hbsc raised level_c into band. At nch=3 ace
% RAISES the slope and hbsc overshoots level_c to 13.9. So this run does NOT
% carry the nch=1 answer over -- it starts from the untouched nch=3 optimum and
% lets the fit find nch=3's own values.
%
% WHAT IS ACTUALLY WRONG WITH nch=3, from that anchor row: it is ALREADY closer
% to the objective than nch=1 ever was (J 0.1152 vs 0.1928), with level_c 4.796
% already IN the [3.5 6.5] band and maperr 147.05 inside its 150 tolerance. The
% whole deficit is two terms:
%     slope 0.5023 -> 0.413    Jslp  0.0893   (78% of J)
%     amp   +37.42 -> >= 40    Jgain 0.0258   (22% of J)
% So maptol stays at the nch=3 default of 150. Unlike the nch=1 run, the ABR here
% is NOT solved, and the map is not the thing to chase -- tightening maptol would
% just pull the search away from the two terms that actually cost something.
%
% 17 PARAMETERS: the default nch=3 fitidx (16, including chsz(3) which nch=1/2 do
% not have, and ace at index 20) plus hbsc via opts.fithbsc. hbsc has never been
% a fitted parameter at ANY chamber count -- it is not in parnames(). gampro
% excluded: CF-map-free but stability-bankrupt at nch=1, and nothing suggests
% more chambers rescue it.
%
% WARM START parfit26_nch3.mat, the best nch=3 fit on record and verified
% sub-critical (maxRe_osc -220.3; all six nch=3 starts on disk were checked
% 2026-08-01 and every one is sub-critical, correcting an earlier claim that this
% one was oscillatory-supercritical at +240).
%
% COST. nch=3 is 184.7 s/eval, 2.6x nch=1. 3 segments x 128 = 384 evals ~ 19.7 h.
% Three segments rather than four because 16-17 parameters need enough
% evaluations per simplex to move at all; two restarts is still enough to tell
% convergence from collapse. That distinction earned its keep on the nch=1 run --
% segment 2 there bought 0.06 of maperr from a FRESH simplex, which is what
% proved segment 1 had genuinely converged rather than stalled.
%
% OUTPUT LOCAL, NOT ONEDRIVE. ~380 checkpoint writes; under OneDrive that is 380
% upload cycles from a conference connection with a sync conflict able to land on
% the live file.
%
% SLEEP. The laptop sleeps in transit Aug 3-9. macOS suspends and resumes MATLAB
% cleanly so the run survives, but the hours printed are WALL clock and read long
% by however long the lid was shut. Judge progress by EVALUATIONS, not hours.
%
% =======================================================================
% ============ WHY THIS PINS d3int, AND WHY IT PRE-FLIGHTS ============
%
% THE FIRST LAUNCH OF THIS SCRIPT WAS INVALID AND BURNED 7.3 h. parfit26 REPORTS
% its start from pa0 but EVALUATES setpar_l(base,pv) with base=modpar26(nch), so
% any tuned field outside the 30 parnames entries + chsz is replaced by a
% modpar26 default on EVERY evaluation while the printed start line still shows
% the good model. The log said `start: slope=0.502 maperr=147.1` and the very
% first evaluation scored map 0.7843, i.e. maperr ~934.
%
% CAUSE (warmdiff.m): parfit26_nch3.mat is from 2026-07-22 and has NO d3int,
% k5*, r5*, m5* fields at all. Today's modpar26(3) ships d3int=1 plus the full
% k5/r5/m5 set. So jointobj was fitting a 3-DOF model while the warm start -- and
% the project's recorded best nch=3 fit, 147.05 -- is 2-DOF. A third DOF switched
% on with parameters never tuned for it is exactly a 147 -> 934 sized error.
%
% pin d3int=0 forces the 2-DOF model onto BOTH base and the warm start, so the
% fit matches its own starting point and the record. k5/r5/m5 go inert with it
% (dof is 3 only when pa.m>=4 || d3int).
%
% The pre-flight below then REFUSES TO RUN if the warm start still does not
% survive setpar_l. This class of defect -- a start report and an evaluation
% disagreeing because they build the model two different ways -- is this
% project's documented dominant failure mode, and it is cheap to assert against.
% =====================================================================
RUNDIR = '/Users/neely/mccm_runs';
if (~exist(RUNDIR,'dir')), mkdir(RUNDIR); end
NSEG = 3; FE = 128;
PIN = struct('d3int',0);

L=load('parfit26_nch3.mat'); pa=L.R.pa;
pf=fieldnames(PIN); for i=1:numel(pf), pa.(pf{i})=PIN.(pf{i}); end

% ---- PRE-FLIGHT: does the warm start survive the parameter mapping? ----
H=parfit26('handles'); b0=modpar26(3);
for i=1:numel(pf), b0.(pf{i})=PIN.(pf{i}); end
pa_eval = H.setpar(b0, H.getpar(pa), []);       % exactly what jointobj builds
try
    e_warm = fdm26(struct('pa',pa));      e_warm = e_warm.maperr;
    e_eval = fdm26(struct('pa',pa_eval)); e_eval = e_eval.maperr;
catch e
    error('c3fit pre-flight could not score: %s', e.message);
end
fprintf('\n  PRE-FLIGHT  warm-start maperr %.2f | as jointobj builds it %.2f | delta %.2f\n', ...
        e_warm, e_eval, abs(e_warm-e_eval));
if (abs(e_warm - e_eval) > 1)
    error(['c3fit ABORTED: the warm start does not survive setpar_l ' ...
           '(maperr %.2f vs %.2f). parfit26 would report one model and fit ' ...
           'another. Pin the missing fields before running.'], e_warm, e_eval);
end
fprintf('  PRE-FLIGHT PASSED -- the model parfit26 fits is the model reported.\n');
fprintf('\n===== C3FIT nch=3: %d segments x %d evals (~%.1f h) =====\n', ...
        NSEG, FE, NSEG*FE*184.7/3600);
fprintf('  warm: parfit26_nch3.mat  ace %+.4f  hbsc %.4f  (J 0.1152)\n', pa.ace, pa.hbsc);
fprintf('  targets: slope 0.5023 -> 0.413 | amp +37.42 -> >=40 | hold lvl_c in band, maperr <=150\n');
fprintf('  17 params: default nch=3 fitidx + hbsc (never fitted at any chamber count)\n');
fprintf('  output: %s\n\n', RUNDIR);

hist=struct('seg',{},'J',{},'slope',{},'lvlc',{},'maperr',{},'amp',{},'osc',{},'ace',{},'hbsc',{},'hrs',{});
bestJ=Inf; bestpa=pa; T0=tic;
for s=1:NSEG
    t0=tic; out=fullfile(RUNDIR,sprintf('c3fit_seg%d.mat',s));
    fprintf('----- segment %d/%d -----\n', s, NSEG);
    try
        R = parfit26(3, struct('maxfe',FE, 'wgain',0.01, 'fithbsc',true, ...
                               'pin',PIN, 'warm',pa, 'verbterm',true, 'out',out));
    catch e
        fprintf('  SEGMENT %d FAILED: %s\n', s, e.message(1:min(150,end)));
        break;
    end
    pa = R.pa;
    S = score26(pa,'fast',false); m = abr_metric(pa,false);
    Js=abs(m.slope-0.413);
    Jl=0.10*(max(0,3.5-m.level_c)+max(0,m.level_c-6.5));
    Jp=0.001*max(0,S.maperr-150); Jg=0.01*max(0,40-S.amp_gain);
    Jo=0.005*max(0,S.maxRe_osc+40);
    J=Js+Jl+Jp+Jg+Jo; hrs=toc(t0)/3600;
    if (J<bestJ), bestJ=J; bestpa=pa; end
    hist(end+1)=struct('seg',s,'J',J,'slope',m.slope,'lvlc',m.level_c, ...
        'maperr',S.maperr,'amp',S.amp_gain,'osc',S.maxRe_osc, ...
        'ace',pa.ace,'hbsc',pa.hbsc,'hrs',hrs); %#ok<SAGROW>
    fprintf('\n  SEG %d | J %.4f | slope %.4f | lvl_c %.3f | maperr %.2f | amp %+.2f | osc %.1f\n', ...
            s, J, m.slope, m.level_c, S.maperr, S.amp_gain, S.maxRe_osc);
    fprintf('        | ace %+.4f | hbsc %.4f | %.1f h\n', pa.ace, pa.hbsc, hrs);
    fprintf('    terms: slope %.4f  lcb %.4f  map %.4f  gain %.4f  osc %.4f\n\n', Js, Jl, Jp, Jg, Jo);
    save(fullfile(RUNDIR,'c3fit_hist.mat'),'hist','bestJ','bestpa');
end

fprintf('\n===== DONE, %.1f h total =====\n', toc(T0)/3600);
fprintf('  seg |      J | slope  | lvl_c | maperr |   amp  |    osc |    ace |  hbsc  | hrs\n');
for i=1:numel(hist)
    h=hist(i);
    fprintf('  %3d | %6.4f | %6.4f | %5.3f | %6.2f | %+6.2f | %6.1f | %+6.4f | %6.4f | %4.1f\n', ...
            h.seg, h.J, h.slope, h.lvlc, h.maperr, h.amp, h.osc, h.ace, h.hbsc, h.hrs);
end
fprintf('\n  best J %.4f    start 0.1152    (nch=1 reached 0.0065)\n', bestJ);
fprintf('  NOTE: the "levers invert at nch=3" premise this run was built on is WRONG.\n');
fprintf('        See the correction block at the top of this file (lcprobe, 2026-08-03).\n');
fprintf('  saved: %s\n', fullfile(RUNDIR,'c3fit_hist.mat'));
disp('C3FIT_DONE');
