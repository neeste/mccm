% CL CARVED FROM SV: the corrected m=4 construction.
%
% SN's decision. CL is carved from SV in BOTH senses at once, and the two are
% the same decision rather than two:
%   AREA      SV gives up CL's share, total conserved.  [0.95 0.05 0.95 0.05]
%   TOPOLOGY  CL takes the position next to the BM, so d1 spans ST-CL, not
%             ST-SV. Stack: ST | BM(d1) | CL | RL(d3) | SS --d2-- SV.
% That is the nested chain. Carving CL out of the SV pool and placing it on the
% BM IS the insert, so pa.nested and the carved chsz belong together.
%
% WHY THE OLD NORMALIZATION CONFUSED THIS. Holding sum(chsz)=2 meant a new
% chamber took its area FROM the existing ones, which is bisection semantics,
% while the topology was an appendage. Mixed. But it took the area from ALL
% chambers by uniform rescale, so it was not a true bisection either. A true
% bisection takes the area from ONE parent and leaves the others untouched.
% Removing the normalization made the areas an appendage too: self-consistent,
% but wrong, since the cochlear cross-section is fixed and resolving CL
% subdivides it. Carving from SV fixes the area half; nesting fixes the
% topology half.
%
% THE 2x2 separates them. Confounding these is exactly what invalidated the
% first nested test (which carried the inverted sign AND the fixed-sum norm).
%
%   row            topology   SV     CL from
%   appended 1.00  appended   1.00   nowhere (total grows to 2.05)  <- built
%   appended 0.95  appended   0.95   SV, area only
%   nested   1.00  nested     1.00   topology only
%   nested   0.95  nested     0.95   SV, both            <- THE CONSTRUCTION
%
% Areas are set explicitly here rather than in modpar26c4 because clsweep2 is
% still running and would hot-reload the edit. The default moves once it lands.

A_APP = [0.95 0.05 1.00 0.05];   % appendage areas, sum 2.05
A_CRV = [0.95 0.05 0.95 0.05];   % SV-carved,       sum 2.00

cfg = { 'm=3b reference  ', 3, [],    0, 0
        'appended SV=1.00', 4, A_APP, 0, 0
        'appended SV=0.95', 4, A_CRV, 0, 0
        'nested   SV=1.00', 4, A_APP, 1, 0
        'nested   SV=0.95', 4, A_CRV, 1, 0
        'nested   SV=0.95 LEGACYSGN', 4, A_CRV, 1, 1 };

fprintf('\n  m=3b ref: amp d1 +81.15  d2 +84.17 | maxRe +19.3 | range 6.53\n');
fprintf('  SN bar for m=4: clean map AND amp 40-60 dB\n\n');
fprintf('  config                     | amp d1  amp d2  amp d3 | maxRe     | range mono fold | maperr\n');
fprintf('%s\n', repmat('-',1,96));
R = struct('nm',{},'S',{});
for i = 1:size(cfg,1)
    pa = modpar26(cfg{i,2});
    if (~isempty(cfg{i,3})), pa.chsz = cfg{i,3}; end
    pa.nested = cfg{i,4}; pa.clvent = 0; pa.m4legacy = cfg{i,5};
    try
        S = score26(pa, 'fast', false);
    catch e
        fprintf('  %-26s | FAILED: %s\n', cfg{i,1}, e.message); continue
    end
    flag = '';
    if (S.maxRe > 24), flag = ' UNSTABLE'; end
    if (S.amp_gain>=40 && S.amp_gain<=60 && strcmp(S.bf_mono,'ok') && S.maxRe<=24)
        flag = ' <== MEETS THE BAR';
    end
    fprintf('  %-26s | %+6.2f %+7.2f %+7.2f | %+9.1f | %5.2f %-4s %.2f | %6.1f%s\n', ...
        cfg{i,1}, S.amp_gain, S.amp_d2, S.amp_d3, S.maxRe, S.bf_range, ...
        S.bf_mono, S.bf_fold, S.maperr, flag);
    R(end+1).nm = cfg{i,1}; R(end).S = S; %#ok<SAGROW>
    save('carve_sv.mat','R');
end
fprintf(['\n  READ AS A 2x2. Compare rows 2-3 for the area carve alone, rows 2 vs 4\n' ...
         '  for the topology alone, and row 5 for both. If row 5 beats every other\n' ...
         '  m=4 row, the construction was the fault and not the parameters. The\n' ...
         '  LEGACYSGN row confirms the sign fix still carries its weight under the\n' ...
         '  new topology rather than being cancelled by it.\n']);
disp('CARVE_SV_DONE');
