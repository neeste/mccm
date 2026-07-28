% RE-TEST: can ohcgain reach 40-60 dB in the 4-chamber, now that score26
% measures the amplifier at the CONFIG'S OWN ohcgain?
%
% WHY THIS MUST RUN BEFORE ANY STRUCTURAL CHANGE. The c4_reasonable sweep
% reported +2.46 dB identically for ohcgain 1, 2 and 4, and I concluded that
% gain cannot raise the m=4 amplifier. That conclusion was UNSOUND for those
% rows: score26 forced ohcgain=1 as the active reference and ignored the swept
% value, so every row was measured at 1. Only the damping rows were valid.
% score26 now uses the config's own ohcgain, so the question is open again.
%
% GRID: n=701 for screening. The convergence test showed AMPLIFIER GAIN and
% maxRe are converged at n=701 (deltas +0.72/+0.30/+1.37 dB; maxRe to 0.1),
% including for the sharply tuned m=3 at +81 dB. maperr and map cleanliness are
% NOT converged there, so those columns are reported but must not be trusted;
% any candidate is re-confirmed at n=1401.
%
% BASE: m2o x32 (clean map) + r1 x0.5, the configuration that met criterion 1
% (fold 0.01, mono ok, maperr 584.8).

b2 = modpar26(4).m2o; br1 = modpar26(4).r1o;
base = modpar26(4); base.m2o = b2*32; base.r1o = br1*0.5;

OG = [1 2 4 8 16 32];
fprintf('\n  n=701 screening (amp and maxRe are converged at this grid)\n');
fprintf('  ohcgain | amp dB | maxRe     | range mono fold | maperr(*)\n');
fprintf('%s\n', repmat('-',1,70));
R = struct('og',{},'S',{});
best = []; bestamp = -Inf;
for og = OG
    pa = setn(base, 701); pa.ohcgain = og;
    try
        S = score26(pa, 'fast', false);
    catch e
        fprintf('  %7g | FAILED: %s\n', og, e.message); continue
    end
    hit = '';
    if (S.amp_gain>=40 && S.amp_gain<=60), hit = '  <== 40-60 dB'; end
    fprintf('  %7g | %+6.2f | %+9.1f | %5.2f %-4s %.2f | %7.1f%s\n', ...
            og, S.amp_gain, S.maxRe, S.bf_range, S.bf_mono, S.bf_fold, S.maperr, hit);
    R(end+1).og=og; R(end).S=S; %#ok<SAGROW>
    save('c4_gain2.mat','R');
    if (isfinite(S.amp_gain) && S.amp_gain > bestamp), bestamp = S.amp_gain; best = og; end
end
fprintf('  (*) maperr at n=701 is NOT converged; shown for trend only.\n');

if (~isempty(best))
    fprintf('\n  Confirming ohcgain=%g at n=1401 (the trustworthy grid)...\n', best);
    pa = base; pa.ohcgain = best;
    S = score26(pa, 'fast');
    save('c4_gain2_confirm.mat','S');
end
fprintf(['\n  If amp still plateaus near 3 dB across a 32x gain range, the limit is\n' ...
         '  structural (the OHC force reaches the BM only as the reaction of an\n' ...
         '  internal pair on d3) and a topology change is justified. If amp rises\n' ...
         '  with ohcgain, the earlier conclusion was a measurement artifact.\n']);
disp('C4_GAIN2_DONE');
