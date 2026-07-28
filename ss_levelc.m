% WHY DOES SS=0.10 BREAK level_c?
%
% THE FINDING. The SS=0.10 refit pre-flight reported level_c = 36.75 at the
% start point, against a target of 5 and a band of [3.5 6.5]. level_c maps to
% per-dB latency shift as %/dB = 100*(level_c^(1/100)-1), so 5 -> 1.62 %/dB (the
% data) and 36.75 -> 3.67 %/dB. The latency shifts with level more than TWICE as
% fast as it should.
%
% WHY THIS WAS MISSED. score26 fast reports maperr, amp, maxRe, range, fold and
% contrast. It never computes level_c. So the basis on which SS-doubling was
% called an improvement on every metric was an INCOMPLETE metric set: it buys an
% 18.5% better maperr (499.3 -> 406.8) and the only clean fold in the project
% (0.84 -> 0.09), and pays in level dependence. The stock [0.95 0.05 1.00] areas
% are FITTED, and part of what they were evidently fitted FOR is level_c.
%
% WHAT THIS SWEEP ANSWERS
%   1. Is the break GRADUAL or a THRESHOLD? A smooth rise from 0.05 to 0.10
%      means SS area is simply a level-dependence lever and the fit balanced it
%      against the map. A sudden jump means something changes state in between,
%      and the place to look is the mode crossing -- fold goes 0.84 -> 0.09 over
%      this same interval, i.e. the BM and shear modes UNCROSS.
%   2. Is there an intermediate SS that keeps most of the map gain while holding
%      level_c in band? That would be the useful configuration, and neither
%      endpoint is it.
%
% THE MODE-CROSSING HYPOTHESIS is the reason to expect a threshold. Uncrossing
% the modes changes which mode carries the peak at the basal places, and the
% level dependence is measured from how the peak MOVES with level. If the peak
% is riding a different mode at high SS, its level behaviour has no reason to
% resemble the fitted one.
%
% COST abr_metric runs tbabr, the expensive path (~150 s per config in the
% pre-flight). Six configs is ~15-20 min, which is why this is worth doing
% directly rather than inferring it from the fit.

SS = [0.05 0.06 0.07 0.08 0.09 0.10];

fprintf('\n  target level_c 5 (=1.62 %%/dB), band [3.5 6.5]\n');
fprintf('  measured at SS=0.10: level_c 36.75 (=3.67 %%/dB), maperr 406.8\n');
fprintf('  fitted   at SS=0.05: maperr 499.3\n\n');
fprintf('  SS    | level_c | %%/dB  | slope | maperr | fold  | in band\n');
fprintf('%s\n', repmat('-',1,64));
R = struct('ss',{},'m',{},'maperr',{},'fold',{});
for s = SS
    pa = modpar26(3); pa.chsz = [0.95 s 1.00];
    lc = NaN; sl = NaN; pdb = NaN;
    try
        evalc('m = abr_metric(pa, false);');
        lc = m.level_c; sl = m.slope;
        if (isfinite(lc) && lc > 0), pdb = 100*(lc^(1/100)-1); end
    catch e
        fprintf('  %.2f  | abr_metric FAILED: %s\n', s, e.message);
    end
    mp = NaN; fo = NaN;
    try
        S = score26(pa, 'fast', false); mp = S.maperr; fo = S.bf_fold;
    catch
    end
    inb = 'no';
    if (isfinite(lc) && lc>=3.5 && lc<=6.5), inb = 'YES'; end
    fprintf('  %.2f  | %7.2f | %5.2f | %5.3f | %6.1f | %5.2f | %s\n', ...
        s, lc, pdb, sl, mp, fo, inb);
    R(end+1).ss=s; R(end).m=lc; R(end).maperr=mp; R(end).fold=fo; %#ok<SAGROW>
    save('ss_levelc.mat','R');
end
fprintf(['\n  READ: a SMOOTH rise means SS area is a level-dependence lever the\n' ...
         '  fit was balancing against the map, and an intermediate value may hold\n' ...
         '  both. A THRESHOLD means something changes state -- compare where the\n' ...
         '  jump sits against where fold collapses (the modes uncrossing). If the\n' ...
         '  two coincide, the level dependence is being measured on a different\n' ...
         '  mode at high SS and the whole comparison is between two regimes.\n']);
disp('SS_LEVELC_DONE');
