% DRIVE THE 4-CHAMBER TO A "REASONABLE" FIT: clean map + 40-60 dB amplifier
% (SN's stated bar, 2026-07-25). The goal is NOT best-fit -- it is a model whose
% inter-chamber displacements and energy flow can be trusted and then mined.
%
% WHERE WE START. m2o x32 already gives a CLEAN map (4.68 oct, mono ok, fold
% removed) with the amplifier intact -- but the amplifier is only +2.46 dB.
% So the map is close to acceptable and THE AMPLIFIER IS THE PROBLEM.
%
% THE QUANTITATIVE TARGET. Amplifier gain tracks maxRe -- proximity to the Hopf
% bifurcation -- monotonically across every model measured, independent of
% chamber count:
%     m=4 m2ox32  maxRe  +0.0 -> amp  +2.46
%     m=2 native  maxRe  +4.5 -> amp +39.11
%     m=1 champ   maxRe  +4.8 -> amp +50.68
%     m=3 native  maxRe +19.3 -> amp +81.15
%     m=3 c3amp   maxRe +29.1 -> amp +116.31
% 40-60 dB therefore corresponds to maxRe ~ +4.5 to +5, exactly where m=1/m=2
% sit. So the task is: RAISE maxRe FROM +0.0 TO ~+5 WITHOUT BREAKING THE MAP.
% (This is the same lesson as the c3amp fit: gain is a SYSTEM property, not a
% local parameter -- so the levers are gain AND damping, i.e. anything that
% moves the system toward criticality.)
%
% CAUTION on ohcgain alone: the 1->0 sweep was near-linear at ~2.4 dB per unit,
% so reaching 40 dB by ohcgain alone would need ~17x and is unlikely to stay
% stable. Damping reduction is the complementary lever and is swept alongside.

b2 = modpar26(4).m2o; base_r1 = modpar26(4).r1o; base_r5 = modpar26(4).r5o;
mk = @(og, fr1, fr5) setfield(setfield(setfield(setfield( ...
        modpar26(4), 'm2o', b2*32), 'ohcgain', og), 'r1o', base_r1*fr1), 'r5o', base_r5*fr5); %#ok<SFLD>

cfg = { 'ohcgain 1  (reference)   ', mk(1,   1,    1)
        'ohcgain 2                ', mk(2,   1,    1)
        'ohcgain 4                ', mk(4,   1,    1)
        'r1 x0.5                  ', mk(1,   0.5,  1)
        'r1 x0.25                 ', mk(1,   0.25, 1)
        'r5 x0.25                 ', mk(1,   1,    0.25)
        'ohcgain 2 + r1 x0.5      ', mk(2,   0.5,  1)
        'ohcgain 4 + r1 x0.25     ', mk(4,   0.25, 1) };

fprintf('\n  config                     | amp dB | maxRe    | range mono fold | maperr\n');
fprintf('%s\n', repmat('-',1,86));
R = struct('tag',{},'S',{});
for i = 1:size(cfg,1)
    tag = cfg{i,1}; pa = cfg{i,2};
    try
        S = score26(pa, 'fast', false);
    catch e
        fprintf('  %-26s | FAILED: %s\n', tag, e.message); continue
    end
    hit = '';
    if (S.amp_gain>=40 && S.amp_gain<=60 && strcmp(S.bf_mono,'ok')), hit = '  <== MEETS THE BAR'; end
    fprintf('  %-26s | %+6.2f | %+8.1f | %5.2f %-4s %.2f | %7.1f%s\n', ...
            tag, S.amp_gain, S.maxRe, S.bf_range, S.bf_mono, S.bf_fold, S.maperr, hit);
    R(end+1).tag = tag; R(end).S = S; %#ok<SAGROW>
    save('c4_reasonable.mat','R');       % bank after each
end
fprintf(['\n  BAR: amp 40-60 dB AND a clean (mono) map. Track maxRe -- it should\n' ...
         '  rise toward ~+5 as amp rises. If amp stays low while maxRe rises, the\n' ...
         '  maxRe/amp relationship does not transfer to m=4 and the 4-chamber gain\n' ...
         '  is limited by its OC-height coupling rather than by criticality.\n']);
disp('C4_REASONABLE_DONE');
