% WHAT DOES ONE m=4 OBJECTIVE EVALUATION ACTUALLY COST?
%
% WHY. The SS=0.10 refit pre-flight measured 148 s per evaluation and I costed a
% 150-evaluation fit at ~6 hours on that basis -- but that was m=3. The
% m4_levelc run implies roughly 30 MINUTES per m=4 abr_metric, an order of
% magnitude more, which would make the same budget ~75 hours. That is the
% difference between a long afternoon and a week, so it needs measuring rather
% than inferring from job timings.
%
% THE BREAKDOWN MATTERS because the two available savings act on different
% parts:
%   skipabr  removes abr_metric (slope + level_c) and keeps fdm26 + coupeig
%   n=701    shrinks everything, but must be shown FAITHFUL, not just faster
% So each part is timed separately at each n, and the METRICS are printed
% alongside: a 2x speedup that moves maperr is not a saving, it is a different
% model. The n=701 convergence test was run earlier in this project; this
% re-checks it specifically for the nested m=4 build, which did not exist then.
%
% CAUTION ON ABSOLUTE NUMBERS. Other jobs may be running, which inflates wall
% times. The RATIOS measured back-to-back inside this script (abr vs map-only,
% 1401 vs 701) are more robust than the absolute values, and the m=3 reference
% row lets the m=4/m=3 ratio be read even under load.
%
% pa.n is NEVER set directly -- setn() resamples the n-length fields (gampro,
% synpro) and rescales isv. Setting pa.n by hand throws or silently corrupts.

NN = [1401 701];

fprintf('\n  Timing one evaluation''s parts. skipabr keeps only the map+stability path.\n');
fprintf('  Metrics are printed so a speedup that MOVES them is visible as such.\n\n');
fprintf('  config          |    n | map+stab s | abr s   | total s | maperr | slope\n');
fprintf('%s\n', repmat('-',1,80));

% ---- m=3 reference, n=1401 only (the 148 s/eval figure came from here) ----
pa = modpar26(3);
t1 = tic; mp = NaN; sl = NaN;
try
    evalc('Rf = fdm26(struct(''pa'',pa));'); mp = Rf.maperr;
    evalc('E = tdm26(''coupeig'', struct(''pa'',pa));'); %#ok<NASGU>
catch e, fprintf('  m=3 map+stab FAILED: %s\n', e.message); end
ta = toc(t1);
t2 = tic;
try, evalc('m = abr_metric(pa,false);'); sl = m.slope; catch, end
tb = toc(t2);
fprintf('  %-15s | %4d | %10.1f | %7.1f | %7.1f | %6.1f | %5.3f\n', ...
    'm=3b reference', pa.n, ta, tb, ta+tb, mp, sl);

% ---- m=4 current default at each n ----
for nn = NN
    pa = modpar26(4);
    if (nn ~= pa.n), pa = setn(pa, nn); end
    mp = NaN; sl = NaN;
    t1 = tic;
    try
        evalc('Rf = fdm26(struct(''pa'',pa));'); mp = Rf.maperr;
        evalc('E = tdm26(''coupeig'', struct(''pa'',pa));'); %#ok<NASGU>
    catch e
        fprintf('  m=4 n=%d map+stab FAILED: %s\n', nn, e.message);
    end
    ta = toc(t1);
    t2 = tic;
    try, evalc('m = abr_metric(pa,false);'); sl = m.slope; catch e2, fprintf('    abr FAILED: %s\n', e2.message); end
    tb = toc(t2);
    fprintf('  %-15s | %4d | %10.1f | %7.1f | %7.1f | %6.1f | %5.3f\n', ...
        'm=4 default', nn, ta, tb, ta+tb, mp, sl);
end

fprintf(['\n  BUDGET ARITHMETIC. Multiply the relevant column by the evaluation\n' ...
         '  count: stage 2a (vent only, 2 params) needs ~40, a full impedance\n' ...
         '  stage ~150. Read the map+stab column for a skipabr fit and the total\n' ...
         '  column for a full one.\n' ...
         '  n=701 IS ONLY USABLE IF maperr AND slope MATCH the n=1401 row. A\n' ...
         '  faster number that differs is a different model, not a saving.\n']);
disp('EVAL_COST_DONE');
