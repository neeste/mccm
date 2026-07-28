% THE RESONANT VENT, NOW WITH maperr -- and the convergence gate vent_res owed.
%
% TWO JOBS IN ONE RUN.
%
% (1) CONVERGENCE GATE. vent_res row B put the CL resonance 4 oct below CF and
% came back 1.5 dB (d1) and 2.6 dB (d2) off row A's eliminated-vent reference.
% At 4 oct below, the vent impedance differs from pure inertance by only ~0.4%
% at CF, so a couple of dB is PLAUSIBLY near-critical sensitivity -- a 2.44%
% area change produced a 54% displacement difference earlier today -- but
% plausible is not verified, and one point cannot separate genuine sensitivity
% from a miswired dof-4 state. Rows B and C push the resonance to 10 and 8 oct
% below CF, where the perturbation is order 1e-6 and 1e-5. If they converge onto
% row A the vent state is sound; if they sit at the same ~2 dB offset as the
% 4-oct row, something is wired wrong and the whole resonance table is void.
%
% (2) maperr FOR THE RESONANT CONFIGS. vent_res measured contrast but not the
% map, because fdm26 had no nested/clvent at the time. The port is now gated and
% passing (regression rows reproduced 1015.2 and 1020.3 exactly), so these are
% the first map measurements of a resonant-vent model.
%
% WHAT vent_res FOUND that this is following up: contrast rises monotonically as
% the CL resonance approaches CF (4.2 at no resonance -> 86.6 at 0.25 oct), but
% the high-contrast rows are unusable (0.50 oct: amp +0.05, maxRe +145.2). The
% usable win was 2.00 oct below CF: contrast 13.2, THREE TIMES the baseline,
% while amp stayed +50.35 (mid-band) and maxRe +10.9 (healthy). The open
% question is what that does to the map, and this answers it.
%
% BASELINE nested + vent 3 -> SS, clvk=0: maperr 522.7, amp +55.08, maxRe +3.9,
% contrast 4.2. m=3b fitted reference: maperr 499.3. nested SEALED: 353.5.

A_CRV = [0.95 0.05 0.95 0.05];
VENT = 3.0; OFF = [10.0 8.0 4.0 2.0 1.0];

base = modpar26(4); base.chsz = A_CRV;
base.nested = 1; base.clvtgt = 2; base.clvent = VENT; base.ohcgain = 1;

try
    evalc('S0 = tdm26(0,base,0,0);');
    cp = S0.cp;
catch e
    fprintf('  could not obtain cp: %s\n', e.message); return
end
w1sq = cp.k1 ./ cp.m1;
clvm = cp.m5 ./ VENT;

fprintf('\n  BASELINE (clvk=0): maperr 522.7 | amp +55.08 | maxRe +3.9 | contrast 4.2\n');
fprintf('  m=3b fitted 499.3 | nested sealed 353.5 | appended 1020.3\n');
fprintf('  Rows B,C are the CONVERGENCE GATE: they must land on row A.\n\n');
fprintf('  config             | maperr  | amp d1  | maxRe    | contrast | range mono\n');
fprintf('%s\n', repmat('-',1,80));
R = struct('lbl',{},'S',{});
for c = 0:numel(OFF)
    pa = base;
    if (c == 0)
        lbl = 'A clvk=0 (ref)';
    else
        pa.clvk = clvm .* (w1sq .* 4^(-OFF(c)));
        lbl = sprintf('%s %.1f oct below', char('A'+c), OFF(c));
    end
    try
        S = score26(pa, 'fast', false);
    catch e
        fprintf('  %-18s | FAILED: %s\n', lbl, e.message); continue
    end
    fprintf('  %-18s | %7.1f | %+6.2f | %+8.1f | %8.1f | %5.2f %-4s\n', ...
        lbl, S.maperr, S.amp_gain, S.maxRe, S.contrast, S.bf_range, S.bf_mono);
    R(end+1).lbl = lbl; R(end).S = S; %#ok<SAGROW>
    save('vent_res2.mat','R');
end
fprintf(['\n  GATE first: do rows B and C sit on row A? Only then read the rest.\n' ...
         '  Then the real question -- does the 2.0 oct row buy its 3x contrast\n' ...
         '  without wrecking the map? If maperr stays near 522.7 there, the\n' ...
         '  resonant vent is a free improvement to tip-tail and the construction\n' ...
         '  meets both halves of SN''s bar with contrast to spare.\n']);
disp('VENT_RES2_DONE');
