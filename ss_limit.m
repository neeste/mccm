% WHAT IS THE STRONG CL->SS VENT LIMIT?
%
% vent_ss saturates: amp d1 +67.27, d2 +80.31, d3 +72.00 at maxRe +3.1, flat
% from vent 30 to vent 100. It does NOT converge to m=3b (d1 is 14 dB short,
% maxRe differs 6x), so my reduction-gate reading was wrong. The identification
% below says what it converges to instead.
%
% THE CLAIM. At p_CL = p_SS the two chambers merge into one and d3 loses its
% driving pressure difference, becoming an internal DOF. That is the m=3b
% TOPOLOGY -- 3 chambers, 3 DOFs, d3 internal -- but with the merged areas
%     ST 0.95,  SS+CL 0.05+0.05 = 0.10,  SV 0.95
% against stock m=3b's 0.95 / 0.05 / 1.00. Different proportions, so the
% numerical mismatch with m=3b is expected rather than a problem.
%
% WHY THIS MATTERS MORE THAN A GATE. If the limit really is a 3-CHAMBER model,
% then fdm26 can score its MAP right now. fdm26 implements neither nested nor
% clvent, which is why maperr has been pinned at the appended value 1020.3 in
% every nested row and why the clean-map half of SN's bar has been untestable.
% The equivalent 3-chamber model has a LIVE maperr. That is a way around the
% blocking item, not merely a consistency check.
%
% ROWS
%   1  m=3b with the merged areas   -> should match the vent-100 row
%   2  m=3b stock areas             -> control, should NOT match (+81.15/+19.3)
%   3  m=3b with SS doubled only    -> separates the SS change from the SV change
%   4  m=3b with SV reduced only    -> the other half of the same separation
% Rows 3 and 4 matter because if row 1 matches, the NEXT question is which of
% the two area changes carries the stability gain (maxRe 19.3 -> 3.1).
%
% TARGET nested vent 100 -> SS: amp d1 +67.27  d2 +80.31  d3 +72.00, maxRe +3.1.

TGT = [67.27 80.31 72.00 3.1];

cfg = { 'm=3b merged  [.95 .10 .95]', [0.95 0.10 0.95]
        'm=3b stock   [.95 .05 1.0]', [0.95 0.05 1.00]
        'm=3b SS only [.95 .10 1.0]', [0.95 0.10 1.00]
        'm=3b SV only [.95 .05 .95]', [0.95 0.05 0.95] };

fprintf('\n  TARGET (nested vent 100 -> SS): d1 %+.2f  d2 %+.2f  d3 %+.2f | maxRe %+.1f\n', TGT);
fprintf('  maperr is LIVE for these rows -- fdm26 supports 3 chambers.\n\n');
fprintf('  config                     | amp d1  amp d2  amp d3 | maxRe     | maperr | range mono\n');
fprintf('%s\n', repmat('-',1,92));
for i = 1:size(cfg,1)
    pa = modpar26(3); pa.chsz = cfg{i,2};
    try
        S = score26(pa, 'fast', false);
    catch e
        fprintf('  %-26s | FAILED: %s\n', cfg{i,1}, e.message); continue
    end
    mk = '';
    if (abs(S.amp_gain-TGT(1))<3 && abs(S.maxRe-TGT(4))<3)
        mk = ' <== MATCHES THE VENT LIMIT';
    end
    fprintf('  %-26s | %+6.2f %+7.2f %+7.2f | %+9.1f | %6.1f | %5.2f %-4s%s\n', ...
        cfg{i,1}, S.amp_gain, S.amp_d2, S.amp_d3, S.maxRe, S.maperr, ...
        S.bf_range, S.bf_mono, mk);
end
fprintf(['\n  IF ROW 1 MATCHES: the CL->SS construction is a 3-chamber model in\n' ...
         '  disguise, its map is scoreable through fdm26 today, and the maperr it\n' ...
         '  reports is the first real measurement of the clean-map half of the bar\n' ...
         '  for any of this session''s configurations. Compare it against m=3b''s\n' ...
         '  499.3 and the appended 1020.3.\n' ...
         '  IF ROW 1 DOES NOT MATCH: the merge identification is wrong, d3 retains\n' ...
         '  some fluid role even at p_CL = p_SS, and nested genuinely does need the\n' ...
         '  fdm26 port before its map can be judged.\n']);
disp('SS_LIMIT_DONE');
