% THE LOW-VENT / HIGH-GAIN CORNER, the one direction vent_cal did not try.
%
% WHAT vent_cal ESTABLISHED. Neither knob alone reaches SN's bar (amp 40-60
% with maxRe <= 24):
%   vent axis at og=1   ratio DEGRADES as vent rises (2.04 -> 1.35 -> 1.15 ->
%                       1.00), so maxRe hits 24 while amp is only ~33.5
%   ohcgain axis at vent 1  a CLIFF between og 0.7 (+11.30, maxRe 6.1) and
%                       og 1.0 (+89.74, maxRe 127.7). Near-critical Hopf.
%
% THE IDEA. Gain per unit instability is best where the vent is WEAKEST (2.04 at
% vent 0.02), and that is also where maxRe has the most headroom: 11.9 against a
% ceiling of 24. So raise gain with ohcgain instead of vent, staying in the
% favourable corner. At ratio 2.04 an amp of 40 costs maxRe ~19.6, inside budget.
%
% THIS EXTRAPOLATION IS OPTIMISTIC AND MAY FAIL. The ratio degrades under
% exactly the kind of pushing proposed here. If it holds near 2.0 the bar is
% reachable; if it collapses toward 1.0 as gain rises, then gain and instability
% are locked together in the nested form regardless of which knob is used, and
% the bar is out of reach without a further structural change.
%
% Rows D re-enter the ohcgain cliff at vent 1.0 between 0.7 and 1.0, to find
% where amp crosses 40-60 and what maxRe costs there. That is the alternative
% route to the bar if the low-vent corner does not work.

A_CRV = [0.95 0.05 0.95 0.05];
SV = 3;

%        label                    vent   ohcgain
cfg = { 'D vent 0.02  og 1.5',    0.02,  1.5
        'D vent 0.02  og 2.0',    0.02,  2.0
        'D vent 0.02  og 3.0',    0.02,  3.0
        'D vent 0.03  og 1.5',    0.03,  1.5
        'D vent 0.03  og 2.0',    0.03,  2.0
        'E vent 1.0   og 0.80',   1.0,   0.80
        'E vent 1.0   og 0.85',   1.0,   0.85
        'E vent 1.0   og 0.90',   1.0,   0.90
        'E vent 1.0   og 0.95',   1.0,   0.95 };

fprintf('\n  BAR: amp d1 in 40-60 AND maxRe <= 24\n');
fprintf('  m=3b 4.20 dB per unit maxRe | nested best so far 2.04 | appended 0.018\n');
fprintf('  headroom at vent 0.02: maxRe 11.9 of 24, amp +24.27\n\n');
fprintf('  config               | amp d1  amp d2  amp d3 | maxRe     | ratio | range mono\n');
fprintf('%s\n', repmat('-',1,86));
R = struct('nm',{},'S',{});
best = []; bestgap = Inf;
for i = 1:size(cfg,1)
    pa = modpar26(4); pa.chsz = A_CRV;
    pa.nested = 1; pa.clvtgt = SV;
    pa.clvent = cfg{i,2}; pa.ohcgain = cfg{i,3};
    try
        S = score26(pa, 'fast', false);
    catch e
        fprintf('  %-20s | FAILED: %s\n', cfg{i,1}, e.message); continue
    end
    rat = NaN;
    if (isfinite(S.amp_gain) && isfinite(S.maxRe) && S.maxRe > 0.1)
        rat = S.amp_gain / S.maxRe;
    end
    flag = '';
    if (S.maxRe > 24), flag = ' unstable'; end
    if (S.amp_gain>=40 && S.amp_gain<=60 && S.maxRe<=24)
        flag = ' <== MEETS THE AMPLIFIER BAR';
        gap = abs(S.amp_gain-50);
        if (gap < bestgap), bestgap = gap; best = cfg{i,1}; end
    end
    fprintf('  %-20s | %+6.2f %+7.2f %+7.2f | %+9.1f | %5.2f | %5.2f %-4s%s\n', ...
        cfg{i,1}, S.amp_gain, S.amp_d2, S.amp_d3, S.maxRe, rat, ...
        S.bf_range, S.bf_mono, flag);
    R(end+1).nm = cfg{i,1}; R(end).S = S; %#ok<SAGROW>
    save('vent_corner.mat','R');
end
if (~isempty(best))
    fprintf('\n  BEST ROW INSIDE THE AMPLIFIER BAR: %s\n', best);
else
    fprintf('\n  NO ROW INSIDE THE AMPLIFIER BAR.\n');
end
fprintf(['  Note the map half of the bar remains UNTESTED for every nested row:\n' ...
         '  fdm26 has no nested or clvent, so maperr cannot respond. Porting that\n' ...
         '  is the blocking item before any nested config can be called good.\n']);
disp('VENT_CORNER_DONE');
