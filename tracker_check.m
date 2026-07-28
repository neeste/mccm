% VERIFY THE CONTINUITY-CONSTRAINED BF TRACKER.
%
% WHAT IT REPLACES. A per-place global argmax, which reported the mode ORDERING
% rather than the map. fold_probe pinned the mechanism exactly: at x/L 0.101 the
% weighted curve peaks at 10.449 kHz with a runner-up at 5.200 kHz just 0.29 dB
% down; one place basal, at x/L 0.094, the same pair reads 5.835 and 10.791 kHz
% separated by 0.01 dB. The peaks swap rank and the argmax reports a 0.84 oct
% cliff across smooth mechanics. The tracker follows the peak nearest in log
% frequency to the previous place instead.
%
% EXPECTED
%   fold      COLLAPSES. The 0.84 (m=3b) and 1.28/1.43 (m=4) folds were mode
%             swaps. If they survive, the tracker is not doing its job.
%   mono      should read ok where fold drops below the 0.10 threshold.
%   range     should GROW: the argmax was pinning to whichever mode won, so the
%             tracked map covers a wider span than the mixture did.
%   maperr    MUST NOT MOVE. It comes from fdm26 and never touched this
%             detector. Movement means the edit reached beyond local_tip.
%   contrast  may shift, since it is read at the tracked peak rather than the
%             global one. A change is legitimate; it is now measured on a
%             consistent mode.
%
% BEFORE (post-edge-fix, pre-tracker):
%   new default  maperr 525.2  range 4.44 FOLD  fold 1.28  chi 6.4  bf 0.34-7.40
%   clvoct off   maperr 522.7  range 4.52 FOLD  fold 1.43  chi 4.2  bf 0.32-7.30
%   m=3b         maperr 499.3  range 5.74 FOLD  fold 0.84  chi 9.2  bf 0.20-10.45
%
% ORIGINAL (pre-edge-fix, both faults present) for reference only:
%   new default  range 6.72 FOLD fold 6.72 | m=3b range 6.53 FOLD fold 6.53
% Those were computed over a set containing bogus apical BFs of 24.98 kHz and
% are not a meaningful baseline for anything.

fprintf('\n  BEFORE (edge fix only, global argmax):\n');
fprintf('    new default  maperr 525.2  range 4.44 fold 1.28  chi 6.4\n');
fprintf('    clvoct off   maperr 522.7  range 4.52 fold 1.43  chi 4.2\n');
fprintf('    m=3b         maperr 499.3  range 5.74 fold 0.84  chi 9.2\n');
fprintf('  maperr must NOT move. fold should COLLAPSE, range should GROW.\n\n');

p1 = modpar26(4);
p2 = modpar26(4); p2 = rmfield(p2,'clvoct');
p3 = modpar26(3);
lbl = {'new default   ', 'clvoct off    ', 'm=3b          '};
PA  = {p1, p2, p3};
fprintf('  config          | maperr  | range mono | bf_lo   bf_hi  | fold  | contrast\n');
fprintf('%s\n', repmat('-',1,80));
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
fprintf(['\n  Fold collapsing with maperr unchanged means the residual folds were\n' ...
         '  mode swaps, not map defects, and the maps in this project have been\n' ...
         '  sound all along. A surviving fold is a REAL reversal and needs\n' ...
         '  chasing on its own.\n']);
disp('TRACKER_CHECK_DONE');
