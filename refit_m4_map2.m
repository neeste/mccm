% m=4 MAP FIT, STABILITY ENFORCED. Rerun of refit_m4_map with the stability
% blind spot closed.
%
% WHAT WENT WRONG LAST TIME. cheapstab=1 removed the per-eval coupeig, so the
% objective could not see stability. The fit bought maperr 376.8 -> 330.6 by
% pushing maxRe +3.7 -> +27.6, past the 24 threshold. Result void.
%
% cheapstab=0 ALONE WOULD NOT HAVE STOPPED IT, and this is the part worth
% getting right. The objective has TWO stability terms and neither would have
% fired:
%   wcrit*max(0, maxRe_OSC + 40)      -- maxRe_osc read -2.3 on the void result,
%                                        i.e. "sub-critical". It filters to
%                                        oscillatory modes above 1 kHz and is
%                                        BLIND to static divergence, which is
%                                        what +27.6 was.
%   wcrit*max(0, maxRe - statol)      -- the static guard, but statol DEFAULTS
%                                        TO 100. At maxRe 27.6 this is zero.
% So the run would have paid 7x the cost for the same unusable answer.
%
% THE FIX
%   cheapstab = 0   coupeig back in the loop (~124 s of the ~144 s per eval)
%   statol    = 20  penalise above 20, below score26's 24 threshold for margin
%   wcrit     = 0.05  sized so instability outweighs any map gain: a 50-point
%                     maperr improvement is worth 0.05 in J, while maxRe at 30
%                     now costs 0.05*10 = 0.5, an order of magnitude more.
%
% COST ~144 s per evaluation against ~20 s before, so roughly 6 hours for 150
% evaluations. That is the price of a search that cannot cheat. Still far below
% the ~4 days a full objective with abr_metric would need.
%
% UNCHANGED: clvent pinned at 0.5 (not in parfit26's parameter space; vent_slope
% measured it as dominating the shipped 3.0 on maperr, level_c and amplifier
% position). chsz pinned to preserve the sum=2 bisection construction, which
% parfit26 cannot enforce as a constraint.
%
% REFERENCES  start 376.8 | void run reached 330.6 at maxRe +27.6
%             m=3b fitted 499.3 | m=1 fitted 104.6
% Expect a stable optimum WORSE than 330.6 -- that number was only reachable by
% going unstable. Anything near 350-400 that holds maxRe under 24 is the real
% result.

pa = modpar26(4);
pa.clvent = 0.5;

opts = struct();
opts.wslope    = 0;
opts.wlevel    = 0;
opts.wlcb      = 0;
opts.skipabr   = 1;
opts.cheapstab = 0;      % <== coupeig per eval
opts.statol    = 20;     % <== was 100; the real guard against static divergence
opts.wcrit     = 0.05;   % <== was 0.005; sized to outweigh map gains
opts.maptol    = 0;
opts.fitidx    = [1 2 3 4 5 7 8 9 10 11 13 20 21];
opts.pin       = struct('chsz', pa.chsz, 'clvent', 0.5, 'clvoct', pa.clvoct, ...
                        'nested', 1, 'clvtgt', 2);
opts.warm      = pa;
opts.maxfe     = 150;
opts.out       = 'refit_m4_map2.mat';

fprintf('\n  m=4 MAP FIT with stability enforced (statol=%g, wcrit=%g, cheapstab=0)\n', ...
    opts.statol, opts.wcrit);
fprintf('  start maperr ~376.8 at maxRe +3.7. Void run hit 330.6 at maxRe +27.6.\n');
fprintf('  A stable optimum should be WORSE than 330.6 -- that was bought with\n');
fprintf('  instability. ~144 s per eval, so ~6 h for %d evals.\n\n', opts.maxfe);

t0 = tic; R = parfit26(4, opts); el = toc(t0);

fprintf('\n  --- m=4 MAP FIT (stability enforced) %.0f s ---\n', el);
if (isfield(R,'pa'))
    fprintf('  n / chsz     : %d / [%s]\n', R.pa.n, strtrim(sprintf('%.3f ', R.pa.chsz)));
    fprintf('  clvent after : %.3f\n', R.pa.clvent);
end
if (isfield(R,'Rf') && isfield(R.Rf,'maperr'))
    fprintf('  final maperr : %.1f   (start 376.8 | m=3b 499.3 | m=1 104.6)\n', R.Rf.maperr);
end
try
    evalc('E = tdm26(''coupeig'', struct(''pa'',R.pa));');
    vd = {'UNSTABLE -- result void', 'healthy'};
    fprintf('  maxRe        : %+.1f (%s)   [osc %+.1f]\n', E.maxRe, ...
        vd{1 + (E.maxRe < 24)}, E.maxRe_osc);
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
fprintf(['\n  Read maxRe FIRST and read the FULL value, not osc. On the void run\n' ...
         '  osc said -2.3 (sub-critical) while the true maxRe was +27.6.\n' ...
         '  amp and contrast remain checks, not targets -- the objective still\n' ...
         '  sees only the map and stability.\n']);
disp('REFIT_M4_MAP2_DONE');
