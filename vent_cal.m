% CALIBRATE THE NESTED + SV-VENT MODEL AGAINST SN's BAR: amp 40-60 dB with a
% clean map, which here means maxRe back inside the healthy band (<=24).
%
% WHERE WE ARE. vent_sv confirmed the choked-series account: venting CL to its
% parent SV restores an ST<->SV route and takes the amplifier from -0.01 dB to
% +115 dB. Two things are now true at once:
%   - the TOPOLOGY question is settled, nested + SV vent works
%   - every vent > 0 is UNSTABLE (maxRe 41.8 at vent 0.1, 255 at vent 3)
% and vent 0.1 already lands amp at +41.69, inside the band, while maxRe sits
% at +41.8, outside it. So the remaining fault is loop GAIN, not topology.
%
% TWO KNOBS, TWO TARGETS. vent sets how well the wave propagates; ohcgain sets
% how hard the OHC drives it. Part A walks the vent below 0.1 to see whether
% propagation alone can land inside both bands. Part B holds a vent that
% propagates well and backs ohcgain off, which is the direct attack on the
% uncalibrated loop gain the bisection predicted.
%
% REFERENCE m=3b reaches +81.15 dB at maxRe +19.3. That is the efficiency to
% beat, and nested currently needs maxRe 127.7 for +89.74. If backing off
% ohcgain moves nested toward m=3b's ratio, the construction is sound and only
% mis-tuned. If amp and maxRe fall together at a FIXED ratio, then nested has a
% worse gain-per-instability constant and no amount of tuning fixes it.
%
% CAVEAT the map half of the bar is NOT testable here: fdm26 has no nested or
% clvent, so maperr is pinned at the appended value 1020.3 in every row. Only
% the tdm26 range column carries any map information at all.

A_CRV = [0.95 0.05 0.95 0.05];
SV = 3;

%        label                     vent   ohcgain
cfg = { 'A vent 0.02  og 1  ',     0.02,  1
        'A vent 0.05  og 1  ',     0.05,  1
        'A vent 0.07  og 1  ',     0.07,  1
        'B vent 1.0   og 0.7',     1.0,   0.7
        'B vent 1.0   og 0.5',     1.0,   0.5
        'B vent 1.0   og 0.3',     1.0,   0.3
        'B vent 1.0   og 0.2',     1.0,   0.2
        'B vent 1.0   og 0.1',     1.0,   0.1
        'C vent 3.0   og 0.3',     3.0,   0.3
        'C vent 3.0   og 0.1',     3.0,   0.1 };

fprintf('\n  BAR: amp d1 in 40-60 AND maxRe <= 24\n');
fprintf('  m=3b: +81.15 at maxRe +19.3 (ratio 4.20 dB per unit maxRe)\n');
fprintf('  nested vent 1 og 1: +89.74 at maxRe +127.7 (ratio 0.70)\n\n');
fprintf('  config               | amp d1  amp d2  amp d3 | maxRe     | ratio | range mono\n');
fprintf('%s\n', repmat('-',1,86));
R = struct('nm',{},'S',{});
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
        rat = S.amp_gain / S.maxRe;   % dB of gain bought per unit instability
    end
    flag = '';
    if (S.maxRe > 24), flag = ' unstable'; end
    if (S.amp_gain>=40 && S.amp_gain<=60 && S.maxRe<=24)
        flag = ' <== MEETS THE AMPLIFIER BAR';
    end
    fprintf('  %-20s | %+6.2f %+7.2f %+7.2f | %+9.1f | %5.2f | %5.2f %-4s%s\n', ...
        cfg{i,1}, S.amp_gain, S.amp_d2, S.amp_d3, S.maxRe, rat, ...
        S.bf_range, S.bf_mono, flag);
    R(end+1).nm = cfg{i,1}; R(end).S = S; %#ok<SAGROW>
    save('vent_cal.mat','R');
end
fprintf(['\n  READ THE RATIO COLUMN. m=3b buys 4.20 dB per unit maxRe. If any\n' ...
         '  nested row approaches that, the construction is efficient and was\n' ...
         '  merely over-driven. If the ratio stays near 0.70 no matter how far\n' ...
         '  ohcgain is backed off, then gain and instability are locked together\n' ...
         '  in the nested form and tuning cannot reach the bar.\n']);
disp('VENT_CAL_DONE');
