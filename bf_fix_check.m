% DID THE APICAL BF FIX WORK, AND WHAT DID IT CHANGE?
%
% THE FAULT. tdm26's find_bf took a GLOBAL max over the whole frequency vector
% with no band limit and no validity test, and the sensitivity curve is
% f-weighted (vh carries an explicit *s factor). At apical places the response
% above CF is negligible, so the weighted curve rises monotonically to the top
% bin and the detector returned the band EDGE. Observed on the seven default
% places, apex to base:
%     BF  = 24.98  24.98  0.78  1.61  3.22  6.37  9.45
%     pbf = -39.07 -38.57  0.54  8.85  6.50  7.43  2.75
% The two most apical places both returned 24.98 kHz -- the last bin, identical
% to the digit -- at levels around -39 dB, against true BFs below 0.78 kHz.
% score26's local_tip has the same f-weighting; being banded to 0.15-18 kHz it
% would pin to the 18 kHz limit rather than 24.98 instead of failing outright.
%
% THE FIX. Reject a maximum sitting on either edge of the admissible range: an
% edge maximum means the curve never turned over, so no BF was found and NaN is
% the correct answer. Applied in BOTH detectors.
%
% WHY IT MATTERS BEYOND TIDINESS. A map reading 24.98, 24.98, 0.78, ... has one
% enormous downward step spanning its whole range, which is exactly the
% fold==range signature that appeared on EVERY model this session including
% m=3b, the accepted scaffolding. I previously called that a metric artifact.
% If this fix clears it, it was a detector failure, and bf_range/bf_fold/bf_mono
% were contaminated for every configuration measured today.
%
% WHAT TO CHECK
%   1. Does the printed BF list still show 24.98 at the apical places?
%   2. Does bf_mono stop reading FOLD?
%   3. bf_range WILL move, because dropping bogus BFs changes max/min. That is a
%      CORRECTION, not a regression -- but it means every bf_range quoted today
%      was computed over a contaminated set.
%   4. maperr must NOT move at all. It comes from fdm26 and never touched these
%      detectors, so any change there means the edit reached further than
%      intended and something else broke.
%
% maperr REFERENCES (must be unchanged): new default 525.2, clvoct off 522.7,
% legacy appended 1015.2, m=3b 499.3.

fprintf('\n  BEFORE (all measured pre-fix):\n');
fprintf('    new default      maperr 525.2  range 6.72 FOLD  contrast 6.4\n');
fprintf('    clvoct off       maperr 522.7  range 6.72 FOLD  contrast 4.2\n');
fprintf('    m=3b             maperr 499.3  range 6.53 FOLD\n');
fprintf('  maperr must NOT move. range/mono SHOULD.\n\n');

fprintf('  --- tdm26 BF list, new default (watch for 24.98 at the apex) ---\n');
pa = modpar26(4);
S0 = tdm26(0, pa, 0, 0); %#ok<NASGU>

fprintf('\n  config          | maperr  | range mono | bf_lo   bf_hi  | fold  | contrast\n');
fprintf('%s\n', repmat('-',1,80));
p1 = modpar26(4);
p2 = modpar26(4); p2 = rmfield(p2,'clvoct');
p3 = modpar26(3);
lbl = {'new default   ', 'clvoct off    ', 'm=3b          '};
PA  = {p1, p2, p3};
for i = 1:numel(PA)
    try
        S = score26(PA{i}, 'fast', false);
    catch e
        fprintf('  %-15s | FAILED: %s\n', lbl{i}, e.message); continue
    end
    fprintf('  %-15s | %7.1f | %5.2f %-4s | %6.2f %6.2f | %5.2f | %8.1f\n', ...
        lbl{i}, S.maperr, S.bf_range, S.bf_mono, S.bf_lo, S.bf_hi, ...
        S.bf_fold, S.contrast);
end
fprintf(['\n  maperr unchanged + mono no longer FOLD = the fix worked and the fold\n' ...
         '  flag was a real detector failure, not a metric quirk. maperr moving\n' ...
         '  means the edit reached somewhere it should not have.\n']);
disp('BF_FIX_CHECK_DONE');
