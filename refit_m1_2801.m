% CAN m=1 BE REFITTED AT n=2801 TO RECOVER ITS MAP?
%
% THE HYPOTHESIS UNDER TEST. Every model's maperr WORSENS under refinement:
%   m=1  70.5 / 79.5 / 88.9 / 104.6 / 130.8 / 155.9 / 211.8   (monotonic)
%   m=3  956.8 / 401.8 / 438.0 / 499.3 / 636.7 / 790.4 / 1114.9  (U, min at 701)
%   m=4  762.5 / 557.0 / 525.1 / 525.2 / 584.7 / 663.9 / 841.1   (U, min ~1401)
% A metric that degrades as resolution improves points at the models. The
% mechanism that fits: coarse grids carry NUMERICAL DAMPING, these models are
% near-critical, so artificial damping suppresses the amplifier. Refining
% removes it, the model runs more active than intended, the peaks shift, and
% maperr worsens. On that account the fitted parameters have ABSORBED the
% numerical damping of the grid they were fitted on -- which is SN's point about
% active forces existing to overcome large passive damping, arriving from the
% numerical side.
%
% THE TEST. If the parameters are merely grid-calibrated, a refit at n=2801
% should recover maperr near the ~100 that m=1 reaches at n=1401. If it plateaus
% near 200, the model cannot represent the data once the numerical damping is
% removed, and the amplifier gains reported throughout this project are partly
% numerical.
%
% WHY m=1. Its advantage is grid-ROBUST (better than m=3/m=4 at every n tested),
% m=1==m=2 holds exactly at every n, and it is cheap. The same test at m=4 would
% be days.
%
% CONFIGURATION, and each choice matters:
%   wslope=0, wlevel=0, wlcb=0   -> needm false at parfit26.m:164, so the
%       per-eval abr_metric is skipped. skipabr alone would NOT do this; it only
%       affects the start and end calls. Slope and level_c are excluded
%       deliberately: both proved unreliable today (slope non-monotonic in n and
%       returning negative values; level_c out of band for every configuration
%       including the fitted scaffolding).
%   cheapstab=1                  -> skips the per-eval coupeig, which at n=2801
%       is an eigenproblem of order dof*n = 5602 and would dominate the cost.
%       Stability is therefore NOT enforced during the search and MUST be
%       checked at the end -- a fit that reaches a good maperr by going unstable
%       is not a result.
%   fitidx without 31-33         -> chsz pinned. m=1's [1 1] symmetry is what
%       makes m=1==m=2 exact; fitting it would break the identity being relied on.
%
% HARD RULE OBSERVED: parfit26 prints the START objective before searching
% (line 118), so the run self-checks that it does not begin from an infeasible
% point. A 7-hour fit was lost to that once.

pa = modpar26(1);
pa = setn(pa, 2801);

opts = struct();
opts.wslope    = 0;      % excludes abr_metric from the per-eval path
opts.wlevel    = 0;
opts.wlcb      = 0;
opts.skipabr   = 1;      % also skip it at start/end
opts.cheapstab = 1;      % no per-eval coupeig at order 5602
opts.fitidx    = [1 2 3 4 5 7 8 9 10 11 13 20 21];
opts.pin       = struct('chsz', [1 1]);
opts.warm      = pa;
opts.maxfe     = 150;
opts.out       = 'refit_m1_n2801.mat';

fprintf('\n  m=1 at n=2801, PURE MAP FIT (slope and level_c excluded).\n');
fprintf('  start maperr should be ~155.9 (n=2801 default); n=1401 default is 104.6.\n');
fprintf('  RECOVERS to ~100 -> parameters were merely grid-calibrated.\n');
fprintf('  PLATEAUS near 200 -> the model cannot fit once numerical damping is gone.\n');
fprintf('  fitting %d params, chsz pinned to [1 1]\n\n', numel(opts.fitidx));

t0 = tic;
R = parfit26(1, opts);
el = toc(t0);

fprintf('\n  --- REFIT SUMMARY (%.0f s) ---\n', el);
if (isfield(R,'Rf') && isfield(R.Rf,'maperr'))
    fprintf('  final maperr : %.1f\n', R.Rf.maperr);
end
if (isfield(R,'pa'))
    fprintf('  chsz after   : [%s]  <== MUST still be 1 1\n', ...
        strtrim(sprintf('%.3f ', R.pa.chsz)));
    fprintf('  n after      : %d\n', R.pa.n);
end
% stability was NOT enforced during the search -- check it now
try
    evalc('E = tdm26(''coupeig'', struct(''pa'',R.pa));');
    vd = {'UNSTABLE -- result void', 'healthy'};
    fprintf('  maxRe        : %+.1f (%s)\n', E.maxRe, vd{1 + (E.maxRe < 24)});
catch e
    fprintf('  coupeig check FAILED: %s\n', e.message);
end
fprintf(['\n  Then re-run the n sweep on the refitted parameters: if maperr is flat\n' ...
         '  across n AFTER refitting, the models are sound but grid-calibrated and\n' ...
         '  every fit simply has to state its grid. If it is still sloped, the\n' ...
         '  refit only moved the minimum and the problem is structural.\n']);
disp('REFIT_M1_2801_DONE');
