% VENT CL TO SS. SN's suggestion, and anatomically the strongest of the three:
% the cortilymph spaces open medially toward the inner sulcus, whereas SV lies
% behind Reissner's membrane and the whole of scala media.
%
% PREDICTION, RECORDED BEFORE THE RUN so it can be wrong. I expect this NOT to
% restore the amplifier, for a structural reason. In the nested chain
%   ST -d1- CL -d3- SS -d2- SV
% CL and SS are ALREADY coupled, by the d3 partition: entries 8 and 14 carry
% -mu3. A CL->SS vent therefore lands in PARALLEL with d3 and drives
% p_CL -> p_SS, collapsing the very pressure difference that drives d3. That is
% the same structural error the ST vent made against d1, moved one link over.
% And it cannot fix the choke: merging CL with SS gives a pool of 0.05+0.05 =
% 0.10 against ST at 0.95, so d1 still faces a chamber ~10x too small and the
% travelling wave still has no ST<->SV route.
%
% WHAT IT MIGHT BE GOOD FOR INSTEAD. Killing d3's fluid drive is exactly what
% strips CL of its independent compartment, so a strong CL->SS vent should
% collapse m=4 toward m=3b. That would give the nested form the REDUCTION GATE
% it currently lacks, since clcouple is inert under nested (that branch never
% reads cc). The m=3b row is included to test convergence, not as a target.
%
% IF I AM WRONG and the amplifier does recover here, then the choked-series
% story is wrong and what matters is simply giving CL somewhere to move.

A_CRV = [0.95 0.05 0.95 0.05];
SS = 2;

cfg = { 'm=3b (collapse target)', 3, 0, 0,    SS
        'appended (reference)  ', 4, 0, 0,    SS
        'nested sealed         ', 4, 1, 0,    SS
        'nested vent 0.5 -> SS ', 4, 1, 0.5,  SS
        'nested vent 1   -> SS ', 4, 1, 1,    SS
        'nested vent 3   -> SS ', 4, 1, 3,    SS
        'nested vent 10  -> SS ', 4, 1, 10,   SS
        'nested vent 30  -> SS ', 4, 1, 30,   SS
        'nested vent 100 -> SS ', 4, 1, 100,  SS };

fprintf('\n  m=3b: amp d1 +81.15  d2 +84.17 | maxRe +19.3 | maperr 499.3\n');
fprintf('  appended SV=0.95: +80.62 | maxRe +4524.8 | nested sealed: -0.01 | +0.0\n\n');
fprintf('  config                 | amp d1  amp d2  amp d3 | maxRe     | range mono | maperr\n');
fprintf('%s\n', repmat('-',1,89));
R = struct('nm',{},'S',{});
for i = 1:size(cfg,1)
    pa = modpar26(cfg{i,2});
    if (cfg{i,2} >= 4), pa.chsz = A_CRV; end
    pa.nested = cfg{i,3}; pa.clvent = cfg{i,4}; pa.clvtgt = cfg{i,5};
    try
        S = score26(pa, 'fast', false);
    catch e
        fprintf('  %-22s | FAILED: %s\n', cfg{i,1}, e.message); continue
    end
    flag = '';
    if (S.maxRe > 24), flag = ' UNSTABLE'; end
    if (S.amp_gain>=40 && S.amp_gain<=60 && strcmp(S.bf_mono,'ok') && S.maxRe<=24)
        flag = ' <== MEETS THE BAR';
    end
    fprintf('  %-22s | %+6.2f %+7.2f %+7.2f | %+9.1f | %5.2f %-4s | %6.1f%s\n', ...
        cfg{i,1}, S.amp_gain, S.amp_d2, S.amp_d3, S.maxRe, S.bf_range, ...
        S.bf_mono, S.maperr, flag);
    R(end+1).nm = cfg{i,1}; R(end).S = S; %#ok<SAGROW>
    save('vent_ss.mat','R');
end
fprintf(['\n  READ TWO WAYS. (1) Does amp d1 leave zero? If yes the choked-series\n' ...
         '  account is wrong. (2) Do the strong-vent rows converge toward the m=3b\n' ...
         '  row? If yes, this vent is a working reduction gate for the nested form\n' ...
         '  even though it is not an amplifier fix, which is worth having on its\n' ...
         '  own since nested currently has no gate at all.\n']);
disp('VENT_SS_DONE');
