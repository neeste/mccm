% DUAL BISECTION OF SV: sum(chsz)=2 BY CONSTRUCTION.
%
% SN's construction, stated precisely: BOTH new chambers are carved from SV, so
% m=4 is SV bisected twice between m=2 and m=4.
%
%   m=2    ST 1.00                        SV 1.00
%   cut 1  ST 1.00  SS s                  SV 1-s
%   cut 2  ST 1.00  SS s   SV 1-s-c   CL c
%   sum  = 1 + s + (1-s-c) + c = 2        IDENTICALLY, for any s and c
%
% The total is conserved by construction rather than by a check or a rescale.
% This also cuts the free area parameters from four to two (s and c), with ST
% pinned at 1 and SV derived -- consistent with the fewer-parameters principle.
%
% WHAT I HAD WRONG. The current default is [0.95 0.05 0.95 0.05]: the CL carve
% is right (SV 1.00 -> 0.95) but SS's 0.05 came from ST, which sits on the FAR
% SIDE of the basilar membrane and cannot be part of bisecting SV. That was
% inherited from modpar26c3's fitted [0.95 0.05 1.00] and I carried it forward
% without questioning it.
%
% RELATED MEASUREMENT, and it is not encouraging for the areas alone. The m=3
% analogue of cut 1, [1.00 0.05 0.95], scored maperr 504.8 against the fitted
% stock [0.95 0.05 1.00] at 499.3 -- a 1% wash, slightly worse. So the area
% BOOKKEEPING was not what the m=3 fit was using that freedom for. What the SS
% sweep did show is that SS SIZE trades CF map against forward-latency slope:
% maperr 499.3 -> 406.8 as SS goes 0.05 -> 0.10, while slope falls 0.400 ->
% 0.185 away from its 0.413 target.
%
% ROWS 2-4 vary the two carve amounts independently, which the old
% parameterisation could not do cleanly because ST was absorbing part of it.
%
% REFERENCE current default [0.95 0.05 0.95 0.05]: maperr 525.2, amp +56.59,
% maxRe +2.3, range 4.44, fold 1.28, chi 6.4.

%        label                     ST    SS    SV    CL
cfg = { 'current  [.95 .05 .95 .05]', 0.95, 0.05, 0.95, 0.05
        'bisect   [1.0 .05 .90 .05]', 1.00, 0.05, 0.90, 0.05
        'bisect   [1.0 .05 .85 .10]', 1.00, 0.05, 0.85, 0.10
        'bisect   [1.0 .10 .85 .05]', 1.00, 0.10, 0.85, 0.05 };

fprintf('\n  REFERENCE current default: maperr 525.2  amp +56.59  maxRe +2.3  chi 6.4\n');
fprintf('  m=3b fitted scaffolding:   maperr 499.3\n');
fprintf('  Rows 2-4 are ST=1 with SS and CL both carved from SV.\n\n');
fprintf('  config                     | sum  | maperr  | amp d1  | maxRe    | range | fold  | chi\n');
fprintf('%s\n', repmat('-',1,92));
R = struct('nm',{},'chsz',{},'S',{});
for i = 1:size(cfg,1)
    cz = [cfg{i,2} cfg{i,3} cfg{i,4} cfg{i,5}];
    if (abs(sum(cz)-2) > 1e-12)
        fprintf('  %-26s | SUM %.4f ~= 2 -- SKIPPED\n', cfg{i,1}, sum(cz)); continue
    end
    pa = modpar26(4); pa.chsz = cz;
    try
        S = score26(pa, 'fast', false);
    catch e
        fprintf('  %-26s | FAILED: %s\n', cfg{i,1}, e.message); continue
    end
    fprintf('  %-26s | %4.2f | %7.1f | %+6.2f | %+8.1f | %5.2f | %5.2f | %5.1f\n', ...
        cfg{i,1}, sum(cz), S.maperr, S.amp_gain, S.maxRe, S.bf_range, ...
        S.bf_fold, S.contrast);
    R(end+1).nm=cfg{i,1}; R(end).chsz=cz; R(end).S=S; %#ok<SAGROW>
    save('dual_bisect.mat','R');
end
fprintf(['\n  Row 2 is the canonical dual bisection and the one to compare against\n' ...
         '  row 1. If it holds or improves, the default should move there: it is\n' ...
         '  the construction SN described, it conserves the total identically, and\n' ...
         '  it drops two free parameters. If it costs materially, the areas are\n' ...
         '  again doing something the fit wants that the geometry does not, which\n' ...
         '  is what the m=3 analogue already hinted at.\n' ...
         '  NOTE these carry m=3b parameters and are UNFITTED, as is row 1, so the\n' ...
         '  comparison is fair even though none of them is a fitted result.\n']);
disp('DUAL_BISECT_DONE');
