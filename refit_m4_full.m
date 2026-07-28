% m=4 COMBINED FIT: impedance + carve amounts, sum(chsz)=2 enforced.
%
% WHAT IS NEW versus refit_m4_map2 (which reached maperr 336.5 at maxRe +9.2):
%   1. chsz is FITTED, not pinned. dual_bisect measured ~60 maperr points in the
%      carve amounts, against the ~40 the impedance parameters yielded. This is
%      where the remaining improvement was located.
%   2. chszderive=3 derives SV so sum(chsz)=2 holds IDENTICALLY through every
%      search step, enforced inside setpar_l. Verified: with the free parameters
%      moved to [0.9229 0.0552 * 0.0503], SV was derived to 0.9716 and the sum
%      held at 2.0000000000; the same values with SV left alone would have
%      summed to 1.9784, a 2.2% drift.
%   3. Balanced stability weights, now that wstat is separate from wcrit:
%        wcrit=0.005  osc term ~0.19, BELOW the map term (~0.34)
%        wstat=0.02   zero at maxRe +9.2, but 0.35 at +27.6
%        statol=10    ~8x the value of the map gain instability would buy
%      refit_m4_map2 shared one wcrit at 0.05, which made the osc term 1.875 --
%      5.6x the map -- so most of that search went into pushing maxRe_osc down
%      rather than fitting the map. Its 336.5 is therefore not the ceiling.
%
% BOTH DUAL-BISECTION READINGS ARE REACHABLE inside the constrained space, and
% the fit is not told which to prefer: "both cuts in SV" is [1, s, 1-s-c, c],
% "one cut in ST and one in SV" is [1-s, s, 1-c, c]. With only sum=2 enforced
% the search can land anywhere between, so WHERE IT LANDS is evidence about the
% construction rather than an assumption baked into it. SN deferred the ST<SV
% anatomical constraint until after the fit, so nothing here imposes it.
%
% index 33 is left IN fitidx deliberately -- parfit26 should auto-drop it and
% say so. That printed line is the confirmation that chszderive took.
%
% cheapstab=0 is REQUIRED: the static guard reads S.maxRe, which only exists if
% coupeig runs per eval. cheapstab=1 would silence the very guard being added.
% Cost ~144 s/eval, so ~8 h for 200 evals.
%
% REFERENCES (n=1401): start 376.8 | impedance-only fit 336.5 at maxRe +9.2
%                      m=3b fitted 499.3 | m=1 fitted 104.6

pa = modpar26(4);
pa.clvent = 0.5;

opts = struct();
opts.wslope     = 0;
opts.wlevel     = 0;
opts.wlcb       = 0;
opts.skipabr    = 1;
opts.cheapstab  = 0;      % required: static guard needs S.maxRe per eval
opts.wcrit      = 0.005;  % osc term back to ~0.19
opts.wstat      = 0.02;   % static guard, independent of wcrit
opts.statol     = 10;
opts.maptol     = 0;
opts.chszderive = 3;      % SV derived -> sum(chsz)=2 identically
opts.fitidx     = [1 2 3 4 5 7 8 9 10 11 13 20 21 31 32 33 34];
opts.pin        = struct('clvent', 0.5, 'clvoct', pa.clvoct, ...
                         'nested', 1, 'clvtgt', 2);   % chsz NOT pinned
opts.warm       = pa;
opts.maxfe      = 200;
opts.out        = 'refit_m4_full.mat';

fprintf('\n  m=4 COMBINED FIT: impedance + carves, sum(chsz)=2 via chszderive=3\n');
fprintf('  Expect a line saying index 33 was dropped from fitidx -- that confirms\n');
fprintf('  chszderive took. Start maperr ~376.8, J ~0.38 + osc ~0.19.\n');
fprintf('  Beat 336.5 (impedance only) and hold maxRe under 24. ~8 h.\n\n');

t0 = tic; R = parfit26(4, opts); el = toc(t0);

fprintf('\n  --- m=4 COMBINED FIT (%.0f s) ---\n', el);
if (isfield(R,'pa'))
    cz = R.pa.chsz;
    fprintf('  chsz         : [%s]\n', strtrim(sprintf('%.4f ', cz)));
    fprintf('  sum(chsz)    : %.10f  (deviation %.2e -- MUST be ~0)\n', sum(cz), abs(sum(cz)-2));
    fprintf('  clvent       : %.3f (pinned 0.500)\n', R.pa.clvent);
    % which dual-bisection reading did it choose? (not imposed -- read off)
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
         '  the map and stability, so a good map with a wrecked amplifier is not\n' ...
         '  a usable model. Where chsz landed is evidence about which bisection\n' ...
         '  reading the data prefers, and was not imposed.\n']);
disp('REFIT_M4_FULL_DONE');
