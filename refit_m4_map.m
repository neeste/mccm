% m=4 PURE MAP FIT at clvent=0.5, n=1401.
%
% THIS IS THE FIT THAT WAS OUT OF REACH ALL SESSION. The blocker was abr_metric
% at ~2290 s per m=4 evaluation, which put 150 evaluations at four days. A pure
% MAP fit avoids it entirely:
%   wslope/wlevel/wlcb = 0 -> needm false at parfit26.m:164, no per-eval
%                             abr_metric (skipabr alone only affects start/end)
%   cheapstab = 1          -> no per-eval coupeig (order dof*n = 4203 here)
%   what remains            -> fdm26 alone, ~20 s. About an hour for 150 evals.
% The m=1 refit at n=2801 just validated this exact configuration: it recovered
% maperr 155.9 -> 98.7 with maxRe +5.0, so the machinery works and the objective
% responds.
%
% maptol = 0 IS ESSENTIAL. The default 167 hinges the map term out of existence
% wherever maperr is already acceptable; the first m=1 refit attempt ran 153
% evaluations against J = 0.0000 for exactly that reason. With maptol=0, J
% tracks maperr directly. Start J should be ~0.377 (= wmap * 376.8), NOT zero.
%
% clvent = 0.5 rather than the shipped 3.0. vent_slope measured 0.5 as dominating
% 3.0 on every column that matters: maperr 376.8 vs 525.2, level_c 5.69 (IN BAND,
% the only configuration all session to manage it) vs 14.60, amp +46.66 vs +56.59
% (both in the 40-60 band), maxRe +3.7 vs +2.3. clvent is NOT in parfit26's
% parameter space -- getpar_l covers impedance 1-30 and chsz 31+ only -- so it is
% PINNED here rather than fitted.
%
% chsz IS PINNED, deliberately. m=4 exposes four area indices (31-34) and
% dual_bisect showed the carve amounts move maperr by ~60 points, so they are
% live parameters. But freeing them lets the search leave the sum(chsz)=2
% bisection construction SN specified, and parfit26 has no way to enforce that
% constraint -- it would need SV derived from the others, which is a change to
% the parameter mapping rather than an option. Fitting the areas is a separate
% follow-up once that derivation exists.
%
% REFERENCE POINTS (all n=1401 unless noted)
%   m=4 clvent=0.5, unfitted   376.8   <- the start
%   m=4 clvent=3.0, unfitted   525.2   (the shipped default)
%   m=3b fitted                499.3
%   m=1  fitted                104.6   (98.7 refitted at n=2801)
% m=1 remains far ahead. The question is how much of m=4's gap is simply that it
% has never been fitted.

pa = modpar26(4);
pa.clvent = 0.5;

opts = struct();
opts.wslope    = 0;
opts.wlevel    = 0;
opts.wlcb      = 0;
opts.skipabr   = 1;
opts.cheapstab = 1;
opts.maptol    = 0;                       % J tracks maperr directly
opts.fitidx    = [1 2 3 4 5 7 8 9 10 11 13 20 21];   % impedance only, no chsz
opts.pin       = struct('chsz', pa.chsz, 'clvent', 0.5, 'clvoct', pa.clvoct, ...
                        'nested', 1, 'clvtgt', 2);
opts.warm      = pa;
opts.maxfe     = 150;
opts.out       = 'refit_m4_map.mat';

fprintf('\n  m=4 PURE MAP FIT, n=%d, clvent=0.5, chsz pinned [%s]\n', ...
    pa.n, strtrim(sprintf('%.2f ', pa.chsz)));
fprintf('  START should report maperr ~376.8 and J ~0.377.\n');
fprintf('  J = 0.0000 at the start means maptol did not take -- STOP, the run is void.\n');
fprintf('  fitting %d impedance params\n\n', numel(opts.fitidx));

t0 = tic; R = parfit26(4, opts); el = toc(t0);

fprintf('\n  --- m=4 MAP FIT SUMMARY (%.0f s) ---\n', el);
if (isfield(R,'pa'))
    fprintf('  n / chsz     : %d / [%s]\n', R.pa.n, strtrim(sprintf('%.3f ', R.pa.chsz)));
    if (isfield(R.pa,'clvent')), fprintf('  clvent after : %.3f (must be 0.500)\n', R.pa.clvent); end
end
if (isfield(R,'Rf') && isfield(R.Rf,'maperr'))
    fprintf('  final maperr : %.1f   (start 376.8 | m=3b 499.3 | m=1 104.6)\n', R.Rf.maperr);
end
% stability was NOT enforced during the search -- verify now
try
    evalc('E = tdm26(''coupeig'', struct(''pa'',R.pa));');
    vd = {'UNSTABLE -- result void', 'healthy'};
    fprintf('  maxRe        : %+.1f (%s)\n', E.maxRe, vd{1 + (E.maxRe < 24)});
catch e
    fprintf('  coupeig check FAILED: %s\n', e.message);
end
% the amplifier was not in the objective either -- check it survived
try
    S = score26(R.pa, 'fast', false);
    fprintf('  amp d1       : %+.2f dB (SN bar 40-60)\n', S.amp_gain);
    fprintf('  contrast     : %.1f | range %.2f\n', S.contrast, S.bf_range);
catch e
    fprintf('  score26 check FAILED: %s\n', e.message);
end
fprintf(['\n  A map-only fit optimises ONE metric. amp and maxRe above are checks,\n' ...
         '  not targets -- if the fit reached a good map by wrecking the amplifier\n' ...
         '  it has not produced a usable model, and that is a real risk given the\n' ...
         '  amplifier was invisible to the objective.\n']);
disp('REFIT_M4_MAP_DONE');
