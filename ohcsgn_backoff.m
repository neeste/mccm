% BACK OFF THE GAIN UNTIL STABLE, with the energy-injecting sign.
%
% Runs only if ohcsgn=-1 was shown to INJECT (ohcP > 0). An injecting pair that
% is unstable is amplifying into self-oscillation, so the remedy is less gain,
% not a different sign.
%
% WHAT "STABLE" MEANS HERE. The champion baseline puts HEALTHY models at
% maxRe = +0.0 to +23.6, while genuinely divergent ones sit at +62747 and above.
% So the target band is maxRe <= ~+24, and anything in the hundreds or thousands
% is self-oscillation regardless of how large the apparent gain is. amp values
% of +219 dB seen earlier were divergence, not amplification: the map had
% collapsed to 0.00 octaves.
%
% GOAL (SN's bar): clean map AND 40-60 dB amplifier. Report the largest gain
% that stays in the healthy maxRe band, and what amplifier it delivers there.
% n=1401 throughout: amp and maxRe are converged at n=701, but map cleanliness
% is NOT, and the map is criterion 1.

b2 = modpar26(4).m2o; br1 = modpar26(4).r1o;
base = modpar26(4); base.m2o = b2*32; base.r1o = br1*0.5; base.ohcsgn = -1;

% sgn=-1 already diverges at ohcgain=1 (maxRe +99507), so the usable range
% is well below 1. Extended down three further decades.
OG = [1 0.5 0.2 0.1 0.05 0.02 0.01 0.005 0.002 0.001];
fprintf('\n  m=4, ohcsgn=-1, m2o x32 + r1 x0.5, n=1401\n');
fprintf('  ohcgain | amp dB | maxRe     | range mono fold | contrast | verdict\n');
fprintf('%s\n', repmat('-',1,80));
R = struct('og',{},'S',{});
bestog = []; bestamp = -Inf;
for og = OG
    pa = base; pa.ohcgain = og;
    try
        S = score26(pa, 'fast', false);
    catch e
        fprintf('  %7g | FAILED: %s\n', og, e.message); continue
    end
    stable = isfinite(S.maxRe) && S.maxRe <= 24;
    clean  = strcmp(S.bf_mono,'ok');
    v = '';
    if (stable && clean),  v = 'stable + clean map'; end
    if (stable && ~clean), v = 'stable, map folded'; end
    if (~stable),          v = 'UNSTABLE'; end
    if (stable && clean && S.amp_gain>=40 && S.amp_gain<=60), v = 'MEETS THE BAR'; end
    fprintf('  %7g | %+6.2f | %+9.1f | %5.2f %-4s %.2f | %8.1f | %s\n', ...
            og, S.amp_gain, S.maxRe, S.bf_range, S.bf_mono, S.bf_fold, S.contrast, v);
    R(end+1).og=og; R(end).S=S; %#ok<SAGROW>
    save('ohcsgn_backoff.mat','R');
    if (stable && isfinite(S.amp_gain) && S.amp_gain > bestamp)
        bestamp = S.amp_gain; bestog = og;
    end
end
if (~isempty(bestog))
    fprintf('\n  Largest STABLE gain: ohcgain=%g giving amp %+.2f dB\n', bestog, bestamp);
else
    fprintf('\n  No configuration in this range was stable (maxRe <= 24).\n');
end
fprintf(['  If the best stable amplifier is far below 40 dB, backing off gain\n' ...
         '  cannot reach the bar and the limit is structural after all.\n']);
disp('OHCSGN_BACKOFF_DONE');
