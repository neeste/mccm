% RE-TEST THE NESTED TOPOLOGY, now that both faults present at the first test
% have been fixed.
%
% WHY THE EARLIER NEGATIVE DOES NOT STAND. The nested chain was tested on
% 2026-07-25 and reported as failing: max|d1| fell 9.46e-08 -> 2.92e-08 even
% sealed, with the apical place breaking. That run carried BOTH faults since
% identified:
%   (1) the BM force sign was INVERTED (+act on the BM instead of -act). The
%       bisection showed correcting it moves the amplifier +3.43 -> ~+78 dB.
%   (2) chsz was renormalized to a FIXED SUM of 2, so m=4's ST and SV were
%       shrunk 2.44% relative to m=3, which near-critical hypersensitivity
%       turned into a 54% displacement difference.
% Both are now fixed, so the nested topology has never been fairly tested.
%
% TOPOLOGIES
%   appended (built)  ST --d1-- SV --d2-- SS --d3-- CL
%   nested            ST --d1-- CL --d3-- SS --d2-- SV
% Only d1 moves between them (ST-SV becomes ST-CL); d2 and d3 are unchanged.
%
% NOTE the nested a2 stamp predates pa.clcouple, so the nested rows have no
% reduction gate. That is a gap to close if nested proves worth pursuing.
%
% ANATOMICAL ASYMMETRY worth keeping in mind when reading this: the spiral
% sulcus is a medial recess, so APPENDING SS to SV is arguably right. Cortilymph
% is unambiguously INTERIOR, between BM and RL, so CL arguably should be
% inserted. The correct construction may be a bisection on top of an appendage
% rather than either uniform scheme.

ISV = [1136 1005 840 655 466 273 80];
cfg = { 'appended, sign OK  ', 0, 0, 0
        'appended, sign LEGACY', 0, 0, 1
        'nested sealed        ', 1, 0, 0
        'nested vent 0.1      ', 1, 0.1, 0
        'nested vent 1        ', 1, 1, 0
        'nested sealed, LEGACY', 1, 0, 1 };

fprintf('\n  m=3b ref: amp d1 +81.15 | maxRe +19.3 | range 6.53\n');
fprintf('  first test (both faults present): appended max|d1| 9.46e-08,\n');
fprintf('  nested sealed 2.92e-08, nested vent1 6.73e-09\n\n');
fprintf('  config                | max|d1|   | amp d1  amp d3 | maxRe     | range mono\n');
fprintf('%s\n', repmat('-',1,78));
for i = 1:size(cfg,1)
    pa = modpar26(4);
    pa.nested = cfg{i,2}; pa.clvent = cfg{i,3}; pa.m4legacy = cfg{i,4};
    pr = pa; pr.isv = ISV;
    md = NaN;
    try
        evalc('S = tdm26(0,pr,0,0);');
        if (all(isfinite(S.d1(:)))), md = max(abs(S.d1(:))); end
    catch
    end
    try
        Q = score26(pa, 'fast', false);
        fprintf('  %-21s | %9.2e | %+6.2f %+7.2f | %+9.1f | %5.2f %-4s\n', ...
            cfg{i,1}, md, Q.amp_gain, Q.amp_d3, Q.maxRe, Q.bf_range, Q.bf_mono);
    catch e
        fprintf('  %-21s | %9.2e | score26 FAILED: %s\n', cfg{i,1}, md, e.message);
    end
end
fprintf(['\n  READ: does nested now hold max|d1| against appended, rather than\n' ...
         '  losing a factor of 3? The LEGACY rows isolate how much of any change\n' ...
         '  is the sign fix versus the topology.\n']);
disp('NESTED_RETEST_DONE');
