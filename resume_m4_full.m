% RESUME the m=4 COMBINED FIT (impedance + carve amounts, sum(chsz)=2 via
% chszderive=3) that was interrupted mid-run.
%
% WHERE IT STOPPED. refit_m4_full.mat.ckpt.mat records iter 72, feval 132 of the
% 200 budget, J=0.5185, saved 28-Jul-2026 09:37:19. The final refit_m4_full.mat
% was never written, so parfit26 did not reach its end reconstruction -- the run
% was cut off, not completed. J had come down from the ~0.567 start (map 0.377 +
% osc 0.19), so the search was live and descending when it died.
%
% WHY A SCRIPT AND NOT A FLAG. parfit26 has no resume path: it always builds pv0
% from the warm start and calls fminsearch fresh. The checkpoint is a crash
% record, not a restart file. fminsearch also does not store its simplex (17
% vertices for 16 free params); only the current best vertex C.x is saved. So an
% exact mid-run continuation is not possible. What IS correct is a WARM RESTART
% from C.x: rebuild the best-so-far pa, hand it back to parfit26 as opts.warm,
% and let fminsearch reinitialise its simplex around that point and keep
% descending. That is the standard way to continue a stalled Nelder-Mead run.
%
% chsz(3) (SV) is re-derived here so sum(chsz)=2 holds at the warm point, exactly
% as chszderive=3 enforces it inside the fit. setpar_l is private to parfit26, so
% its mapping (30 impedance params by name, chsz in 31+, index 3 derived) is
% replicated below.
%
% COST ~144 s/eval, so the maxfe=120 below is ~5 h. Lower it if you only want to
% top up the ~68 evals the original budget had left.
%
% REFERENCES (n=1401): start 376.8 | impedance-only 336.5 | m=3b 499.3 | m=1 104.6

CK = load('refit_m4_full.mat.ckpt.mat'); C = CK.C;
fprintf('\n  RESUME m=4 combined fit from checkpoint\n');
fprintf('  checkpoint : iter %d, feval %d/200, J %.4f, saved %s\n', ...
        C.iter, C.fevals, C.J, C.when);

% --- rebuild the current-best pa from the checkpoint (replicates setpar_l) ---
pa = modpar26(4);
pa.clvent = 0.5;
nm = {'k1o','r1o','m1o','k2o','r2o','m2o','k3o','r3o','k4o','aco', ...
      'k1e','r1e','m1e','k2e','r2e','m2e','k3e','r3e','k4e','ace', ...
      'k1q','r1q','m1q','k2q','r2q','m2q','k3q','r3q','k4q','acq'};
pv = C.pv0; pv(C.fitidx) = C.x;          % apply the fitted params onto pv0
for i = 1:30, pa.(nm{i}) = pv(i); end
nc = numel(pa.chsz);
pa.chsz = pv(31:30+nc);
o = true(1,nc); o(3) = false;            % derive SV (index 3) so sum(chsz)=2
pa.chsz(3) = 2 - sum(pa.chsz(o));
fprintf('  warm chsz  : [%s] sum %.10f\n', ...
        strtrim(sprintf('%.4f ', pa.chsz)), sum(pa.chsz));

% --- identical objective configuration to refit_m4_full.m ---
opts = struct();
opts.wslope     = 0;
opts.wlevel     = 0;
opts.wlcb       = 0;
opts.skipabr    = 1;
opts.cheapstab  = 0;      % required: static guard needs S.maxRe per eval
opts.wcrit      = 0.005;  % osc term ~0.19
opts.wstat      = 0.02;   % static guard, independent of wcrit
opts.statol     = 10;
opts.maptol     = 0;
opts.chszderive = 3;      % SV derived -> sum(chsz)=2 identically
opts.fitidx     = [1 2 3 4 5 7 8 9 10 11 13 20 21 31 32 33 34];
opts.pin        = struct('clvent', 0.5, 'clvoct', pa.clvoct, ...
                         'nested', 1, 'clvtgt', 2);   % chsz NOT pinned
opts.warm       = pa;     % <== warm restart from the checkpoint's best point
opts.maxfe      = 120;    % remaining budget; ~144 s/eval => ~5 h
opts.out        = 'refit_m4_full.mat';

fprintf('  expect a line saying index 33 was dropped from fitidx (chszderive)\n');
fprintf('  start maperr should read ~ the checkpoint level, NOT 376.8, since the\n');
fprintf('  warm start is the partially-fitted point. Beat 336.5, hold maxRe < 24.\n\n');

t0 = tic; R = parfit26(4, opts); el = toc(t0);

fprintf('\n  --- m=4 COMBINED FIT RESUMED (%.0f s) ---\n', el);
if (isfield(R,'pa'))
    cz = R.pa.chsz;
    fprintf('  chsz         : [%s]\n', strtrim(sprintf('%.4f ', cz)));
    fprintf('  sum(chsz)    : %.10f  (deviation %.2e -- MUST be ~0)\n', sum(cz), abs(sum(cz)-2));
    fprintf('  clvent       : %.3f (pinned 0.500)\n', R.pa.clvent);
    if (abs(cz(1)-1) < 0.01)
        rd = 'ST~1.00: reading A, both cuts in SV';
    elseif (abs(cz(1)-cz(3)) < 0.02)
        rd = 'ST~SV: reading B, one cut in each';
    else
        rd = 'neither reading -- the data prefer something between';
    end
    fprintf('  ST vs SV     : %.4f vs %.4f  -> %s\n', cz(1), cz(3), rd);
end
if (isfield(R,'Rf') && isfield(R.Rf,'maperr'))
    fprintf('  final maperr : %.1f   (start 376.8 | impedance-only 336.5 | m=3b 499.3 | m=1 104.6)\n', R.Rf.maperr);
end
try
    evalc('E = tdm26(''coupeig'', struct(''pa'',R.pa));');
    vd = {'UNSTABLE -- result void', 'healthy'};
    fprintf('  maxRe        : %+.1f (%s)  [osc %+.1f]\n', E.maxRe, vd{1+(E.maxRe<24)}, E.maxRe_osc);
catch e
    fprintf('  coupeig FAILED: %s\n', e.message);
end
try
    S = score26(R.pa, 'fast', false);
    fprintf('  amp d1       : %+.2f dB (bar 40-60) | contrast %.1f | range %.2f\n', ...
        S.amp_gain, S.contrast, S.bf_range);
catch e
    fprintf('  score26 FAILED: %s\n', e.message);
end
fprintf(['\n  amp and contrast are CHECKS, not targets -- the objective sees only\n' ...
         '  the map and stability. Where chsz landed is evidence about which\n' ...
         '  bisection reading the data prefer, and was not imposed.\n']);
disp('RESUME_M4_FULL_DONE');
