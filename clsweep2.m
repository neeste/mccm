% WHAT DOES GIVING d3 A COMPARTMENT ACTUALLY DO?  (re-run after the chsz
% normalization change, which moved m=4's areas 2.44% and which near-critical
% hypersensitivity amplifies. The pre-change sweep is stale.)
%
% pa.clcouple in [0,1] scales d3's fluid coupling. The reduction gate at 0 is
% now EXACT to roundoff (d1 rel diff 3.6e-10 vs m=3b), so this sweep starts from
% a verified baseline and every departure is attributable to the compartment.
%
% This measures the SECOND fault the bisection found. The first, the inverted BM
% force sign, is now fixed by construction (the gate forces -act on the BM).
%
% References: m=3b amp d1 +81.15, d2 +84.17, maxRe +19.3, maperr 499.3.
% SN's bar for m=4: clean map AND 40-60 dB amplifier.
% maperr is reported but fdm26 has NO clcouple, so it cannot respond to the
% sweep; it is shown only to confirm it stays put.

CC = [0 0.02 0.05 0.1 0.2 0.35 0.5 0.7 0.85 1.0];
fprintf('\n  m=3b ref: amp d1 +81.15  d2 +84.17 | maxRe +19.3 | range 6.53\n');
fprintf('  SN bar  : clean map AND amp 40-60 dB\n\n');
fprintf('  clcouple | amp d1  amp d2  amp d3 | maxRe     | range mono fold | maperr\n');
fprintf('%s\n', repmat('-',1,80));
R = struct('cc',{},'S',{});
for cc = CC
    pa = modpar26(4); pa.clcouple = cc;
    try
        S = score26(pa, 'fast', false);
    catch e
        fprintf('  %8.2f | FAILED: %s\n', cc, e.message); continue
    end
    flag = '';
    if (S.maxRe > 24),  flag = ' UNSTABLE'; end
    if (S.amp_gain>=40 && S.amp_gain<=60 && strcmp(S.bf_mono,'ok') && S.maxRe<=24)
        flag = ' <== MEETS THE BAR';
    end
    fprintf('  %8.2f | %+6.2f %+7.2f %+7.2f | %+9.1f | %5.2f %-4s %.2f | %6.1f%s\n', ...
        cc, S.amp_gain, S.amp_d2, S.amp_d3, S.maxRe, S.bf_range, S.bf_mono, ...
        S.bf_fold, S.maperr, flag);
    R(end+1).cc=cc; R(end).S=S; %#ok<SAGROW>
    save('clsweep2.mat','R');
end
fprintf(['\n  Read where maxRe leaves the healthy band (<=24). That value of\n' ...
         '  clcouple is how much fluid coupling d3 can carry before the loop gain\n' ...
         '  runs away, measured from a verified baseline for the first time.\n']);
disp('CLSWEEP2_DONE');
