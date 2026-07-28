% REFIT m=1 AT n=2801 -- corrected setup.
%
% THE FIRST ATTEMPT WAS VOID, for two reasons, both mine.
%
% (1) IT RAN AT THE WRONG GRID. parfit26.m:92 builds `base = modpar26(nch)` at
% the DEFAULT n, and the objective at line 113 closes over `base`. opts.warm
% supplies parameter VALUES only -- the grid comes from base. So passing a
% setn'd pa as warm did nothing and the fit ran at n=1401. Reported "n after:
% 1401", and a final maperr of 104.6, which is exactly the n=1401 default.
%
% FIX: opts.pin forces fields onto BOTH base and the warm start
% (parfit26.m:99-103). But n cannot be pinned alone -- setn() exists because
% changing n requires RESAMPLING the n-length fields (gampro, synpro) and
% rescaling isv; setting pa.n by hand corrupts them silently. So setn is used
% first and every field it touches is pinned together.
%
% (2) THE OBJECTIVE HAD NO GRADIENT. maptol defaults to 167 and the map term is
% wmap*max(0, maperr - maptol). At maperr 104.6 -- or 155.9, the n=2801 default
% -- that term is ZERO, so J was 0.0000 from the first iteration and 153
% evaluations minimised a constant. maptol=0 makes J proportional to maperr
% itself, which is what a pure map fit needs.
%
% This is the second time today an objective was mis-specified in a way that
% wasted a run: the SS=0.10 pre-flight found the default weighting drives maperr
% the WRONG way (wmap=0.001 against a level_c term worth ~3.0). The lesson both
% times is the same -- read the objective's value AND its gradient at the start
% point, not just its value.
%
% WHAT THIS TESTS, at the lower stakes SN set: whether the parameterisation has
% enough freedom to reach a good map at a finer grid. Not model validity --
% parameters absorbing grid-dependent numerical damping is a calibration
% property, and that is accepted.
%   START  maperr 155.9 (m=1 default at n=2801)
%   TARGET ~104.6, what m=1 reaches at n=1401
%   RECOVERS -> grid choice is free, state it and move on
%   PLATEAUS -> the reachable map quality depends on the grid

p1401 = modpar26(1);
p2801 = setn(p1401, 2801);

% pin every field setn touches, or n changes without its dependents
pin = struct('n', p2801.n, 'chsz', [1 1]);
for fn = {'gampro','synpro','isv'}
    if (isfield(p2801, fn{1})), pin.(fn{1}) = p2801.(fn{1}); end
end

opts = struct();
opts.wslope    = 0;      % keeps needm false -> no per-eval abr_metric
opts.wlevel    = 0;
opts.wlcb      = 0;
opts.skipabr   = 1;
opts.cheapstab = 1;      % coupeig at order 5602 would dominate; checked at end
opts.maptol    = 0;      % <== THE FIX: make the map term proportional to maperr
opts.fitidx    = [1 2 3 4 5 7 8 9 10 11 13 20 21];
opts.pin       = pin;
opts.warm      = p2801;
opts.maxfe     = 150;
opts.out       = 'refit_m1_n2801b.mat';

fprintf('\n  pinning: %s\n', strjoin(fieldnames(pin)', ', '));
fprintf('  START should report maperr ~155.9 and n=2801. If it reports 104.6 or\n');
fprintf('  n=1401 the grid pin has failed again and the run is void -- stop.\n');
fprintf('  maptol=0 so J tracks maperr directly (J was identically 0 last time).\n\n');

t0 = tic; R = parfit26(1, opts); el = toc(t0);

fprintf('\n  --- REFIT SUMMARY (%.0f s) ---\n', el);
ok = true;
if (isfield(R,'pa'))
    fprintf('  n after      : %d\n', R.pa.n);
    if (R.pa.n ~= 2801), fprintf('  *** GRID PIN FAILED -- result void ***\n'); ok = false; end
    fprintf('  chsz after   : [%s]\n', strtrim(sprintf('%.3f ', R.pa.chsz)));
end
if (isfield(R,'Rf') && isfield(R.Rf,'maperr'))
    fprintf('  final maperr : %.1f   (start 155.9, n=1401 reference 104.6)\n', R.Rf.maperr);
end
if (ok)
    try
        evalc('E = tdm26(''coupeig'', struct(''pa'',R.pa));');
        vd = {'UNSTABLE -- result void', 'healthy'};
        fprintf('  maxRe        : %+.1f (%s)\n', E.maxRe, vd{1 + (E.maxRe < 24)});
    catch e
        fprintf('  coupeig check FAILED: %s\n', e.message);
    end
end
disp('REFIT_M1_2801B_DONE');
