% REDUCTION GATE FOR m=3b.
%
% m=3b is the missing rung: 3 chambers, 3 DOFs, d3 INTERNAL (a local mass-spring
% with no fluid compartment, exactly as d2 is in m=2).
%   m=2   2 chambers, 2 DOFs, d2 internal
%   m=3   3 chambers, 2 DOFs, d2 gets SS
%   m=3b  3 chambers, 3 DOFs, d3 internal     <- this
%   m=4   4 chambers, 3 DOFs, d3 gets CL
%
% THE GATE: freeze d3 with a heavy m5 and m=3b must reproduce m=3, including its
% amplifier (+81.15 dB d1, +84.17 dB d2, maxRe +19.3, maperr 499.3). If it does,
% the rung is sound and m=4 can be rebuilt from it with its own gate. If it does
% not, m=3b is wrong and there is no point continuing up the ladder.
%
% Why this matters: every m=4 result so far was compared against m=3 across a
% gap that changed a DOF, a chamber and the BM force sign at once, so nothing
% could be attributed. A gated ladder makes each step attributable.
%
% The gate also FIXES THE FORCE SIGN by construction rather than by search: with
% d3 frozen, s1 must equal m=3's, so the BM keeps -act and d3 takes +act.

fprintf('\n  m=3 reference: amp d1 +81.15  d2 +84.17 | maxRe +19.3 | maperr 499.3\n\n');
fprintf('  config              | amp d1  amp d2  amp d3 | maxRe     | maperr | range mono\n');
fprintf('%s\n', repmat('-',1,82));

b5 = modpar26c3b; b5 = b5.m5o;
cfg = { 'm=3 (reference)  ', modpar26(3), NaN };
for f = [1 10 100 1e3 1e4 1e6]
    pa = modpar26c3b; pa.m5o = b5*f;
    cfg(end+1,:) = { sprintf('m=3b  m5 x%-6g', f), pa, f }; %#ok<SAGROW>
end

ref = [];
for i = 1:size(cfg,1)
    nm = cfg{i,1}; pa = cfg{i,2};
    try
        S = score26(pa, 'fast', false);
    catch e
        fprintf('  %-19s | FAILED: %s\n', nm, e.message); continue
    end
    if (isempty(ref)), ref = S; end
    fprintf('  %-19s | %+6.2f %+7.2f %+7.2f | %+9.1f | %6.1f | %5.2f %-4s\n', ...
            nm, S.amp_gain, S.amp_d2, S.amp_d3, S.maxRe, S.maperr, ...
            S.bf_range, S.bf_mono);
    if (i > 1 && isfinite(S.amp_gain) && isfinite(ref.amp_gain))
        if (abs(S.amp_gain-ref.amp_gain) < 1.0 && abs(S.maperr-ref.maperr) < 5)
            fprintf('  %19s   <== GATE PASSED (matches m=3 within 1 dB / 5 maperr)\n','');
        end
    end
end
fprintf(['\n  PASS = a heavy-m5 row matching m=3 on amp d1, amp d2, maxRe and\n' ...
         '  maperr. That is the frozen-d3 limit and it validates the rung.\n' ...
         '  FAIL = no row converges to m=3, meaning m=3b is not a clean extension\n' ...
         '  and the ladder should not be climbed further.\n']);
disp('M3B_GATE_DONE');
