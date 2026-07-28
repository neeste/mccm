% SHOULD m=3's AREAS FOLLOW THE BISECTION TOPOLOGY?
%
% THE ARGUMENT. With the nested rebuild, m=4 is now a DUAL BISECTION of m=2,
% both cuts falling in the space above the BM:
%     m=2    ST --d1-- SV'                      (SV' = everything above the BM)
%     cut 1  SV' -> CL + rest, divided by d3     ST --d1-- CL --d3-- (rest)
%     cut 2  rest -> SS + SV, divided by d2      ST --d1-- CL --d3-- SS --d2-- SV
%
% SN's CL carve is consistent with that: SV 1.00 -> 0.95 with CL taking the
% 0.05, cut and area both in the upper chamber. But the m=2 -> m=3 step is NOT:
%     m=2  ST 1.00  SV' 1.00
%     m=3  ST 0.95  SS 0.05  SV 1.00     <- SS's area came from ST
% ST sits on the FAR SIDE of the basilar membrane. A bisection of the upper
% chamber cannot draw area from the lower one. modpar26c3 carries
% chsz=[0.95 0.05 1.0] with the comment err=93.36, so it is a FITTED value, not
% a construction-derived one, and nothing forced it to respect a bisection that
% had not been articulated at the time.
%
% THE TEST. Make the area follow the cut: take SS's 0.05 from SV instead of ST.
% Same total, same chambers, same topology -- only the bookkeeping changes.
%
% WHY EXPECT ANYTHING. ss_limit already measured this model as sensitive on the
% SS axis: [0.95 0.10 1.0] gave maperr 406.8 against the stock 499.3, an 18%
% improvement from SS's area alone. Nobody has tried the bisection-consistent
% assignment, and row 4 combines both.
%
% CAUTION m=3b is the SCAFFOLDING every comparison in this project is anchored
% to, including the fitted maperr 499.3. A better row here does NOT license
% changing modpar26c3 -- these parameters were fitted around the stock areas, so
% a fair comparison needs a refit at the new areas before anything is adopted.
%
% METRICS come from the reverted detector: argmax plus the apical edge guard.
% The continuity tracker was tried and REVERTED (it followed the sub-dominant
% mode past the crossing, to rank 11 and -16 dB). bf_fold measures MODE
% DEGENERACY here, not map quality -- m=3b's 0.84 oct is the BM and shear modes
% crossing at x/L 0.10 with a 0.01 dB margin.

cfg = { 'stock     [0.95 0.05 1.00]', [0.95 0.05 1.00]
        'bisection [1.00 0.05 0.95]', [1.00 0.05 0.95]
        'SS x2 stk [0.95 0.10 1.00]', [0.95 0.10 1.00]
        'SS x2 bis [1.00 0.10 0.90]', [1.00 0.10 0.90] };

fprintf('\n  REFERENCE m=3b stock: maperr 499.3 | amp +81.15 | maxRe +19.3\n');
fprintf('  ss_limit measured [0.95 0.10 1.00] at maperr 406.8 (SS area alone)\n');
fprintf('  All rows keep sum(chsz)=2.00 and the same topology.\n\n');
fprintf('  config                     | maperr  | amp d1  | maxRe    | range | fold  | contrast\n');
fprintf('%s\n', repmat('-',1,90));
best = Inf; bl = '';
for i = 1:size(cfg,1)
    pa = modpar26(3); pa.chsz = cfg{i,2};
    try
        S = score26(pa, 'fast', false);
    catch e
        fprintf('  %-26s | FAILED: %s\n', cfg{i,1}, e.message); continue
    end
    mk = '';
    if (isfinite(S.maperr) && S.maperr < best), best = S.maperr; bl = cfg{i,1}; end
    if (isfinite(S.maperr) && S.maperr < 499.3 - 5), mk = ' <== better than stock'; end
    fprintf('  %-26s | %7.1f | %+6.2f | %+8.1f | %5.2f | %5.2f | %8.1f%s\n', ...
        cfg{i,1}, S.maperr, S.amp_gain, S.maxRe, S.bf_range, S.bf_fold, ...
        S.contrast, mk);
end
fprintf('\n  best maperr %.1f at %s\n', best, strtrim(bl));
fprintf(['  Row 2 beating row 1 would say the areas SHOULD follow the cut, and\n' ...
         '  that the fitted 0.95/1.00 split was compensating for a topology the\n' ...
         '  model did not have until the nested rebuild. Row 2 losing would say\n' ...
         '  the area bookkeeping is not what the fit was using that freedom for.\n' ...
         '  Either way this needs a REFIT at the new areas before adoption --\n' ...
         '  every parameter in modpar26c3 was fitted around the stock split.\n']);
disp('BISECT_AREA_DONE');
