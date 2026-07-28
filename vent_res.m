% RESONANT VENT: reduction gate first, then place the CL resonance.
%
% WHAT WAS ADDED. The vent was already a pure inertance -- its a2 stamp
% G = clvent*mu3 = clvent*m1/m5 is exactly the coefficient m1/mv of a DOF with
% mass mv = m5/clvent, algebraically eliminated for want of stiffness or
% damping. One number therefore set BOTH the coupling strength and the
% resonance, which is why twomap2 found the realized d1->d3 transfer peak pinned
% at CF regardless of k5 (medians -0.060 vs -0.114 oct for partition resonances
% a full octave apart, indistinguishable at IQR ~0.5).
%
% pa.clvk gives the cortilymph space a compliance, making the vent a Helmholtz
% resonator: mv the channel inertance, clvk the compliance, clvr the channel
% resistance. The resonance sqrt(clvk/mv) is now placeable INDEPENDENTLY of the
% coupling strength, which is the control that was missing.
%
% THE GATE. clvk = clvr = 0 keeps dof at 3 and the vent algebraically eliminated
% exactly as before, so row A must reproduce vent_ss's vent-3 row (+55.08 dB,
% maxRe +3.9) EXACTLY. Row B puts the resonance ~4 oct below CF, far below
% anything that matters, and must come back very close to row A -- that is the
% real test, since it exercises the new dof-4 state path and must agree with the
% eliminated form in the limit where the stiffness is negligible.
% If row B disagrees with row A, the DOF is wired wrong and NOTHING below it
% means anything.
%
% clvk IS PLACE-DEPENDENT. To sit a fixed octave interval below the local CF it
% must track w1(x) = sqrt(k1/m1), so it is derived from cp returned by a first
% run rather than from a guessed x grid.

A_CRV = [0.95 0.05 0.95 0.05];
VENT = 3.0; OFF = [4.0 2.0 1.0 0.5 0.25];

base = modpar26(4); base.chsz = A_CRV;
base.nested = 1; base.clvtgt = 2; base.clvent = VENT; base.ohcgain = 1;

fprintf('\n  REFERENCE vent_ss vent 3 -> SS: amp d1 +55.08  d2 +69.78  d3 +59.52, maxRe +3.9\n');
fprintf('  Row A must reproduce it. Row B (4 oct below CF) must nearly reproduce it.\n\n');

% one run to obtain cp on the model grid
try
    evalc('S0 = tdm26(0,base,0,0);');
    cp = S0.cp;
catch e
    fprintf('  could not obtain cp: %s\n', e.message); return
end
w1sq = cp.k1 ./ cp.m1;            % local BM radian resonance, squared
clvm = cp.m5 ./ VENT;             % channel inertance implied by clvent
fprintf('  grid ok: n=%d, w1/2pi base %.0f Hz apex %.0f Hz\n\n', ...
    numel(w1sq), sqrt(w1sq(1))/(2*pi), sqrt(w1sq(end))/(2*pi));

fprintf('  config            | amp d1  amp d2  amp d3 | maxRe    | contrast | range mono\n');
fprintf('%s\n', repmat('-',1,84));
R = struct('lbl',{},'S',{});
for c = 0:numel(OFF)
    pa = base;
    if (c == 0)
        lbl = 'A clvk=0 (gate)';
    else
        off = OFF(c);
        pa.clvk = clvm .* (w1sq .* 4^(-off));   % resonance off octaves below CF
        lbl = sprintf('%s %.2f oct below', char('A'+c), off);
    end
    try
        S = score26(pa, 'fast', false);
    catch e
        fprintf('  %-17s | FAILED: %s\n', lbl, e.message); continue
    end
    fprintf('  %-17s | %+6.2f %+7.2f %+7.2f | %+8.1f | %8.1f | %5.2f %-4s\n', ...
        lbl, S.amp_gain, S.amp_d2, S.amp_d3, S.maxRe, S.contrast, ...
        S.bf_range, S.bf_mono);
    R(end+1).lbl = lbl; R(end).S = S; %#ok<SAGROW>
    save('vent_res.mat','R');
end
fprintf(['\n  READ ROW B FIRST. If it differs materially from row A, the dof-4 vent\n' ...
         '  state is wired wrong and the rest of the table is meaningless. Only if\n' ...
         '  the gate holds does the contrast column test SN''s hypothesis, and the\n' ...
         '  question there is whether contrast rises as the CL resonance approaches\n' ...
         '  0.5 oct below CF -- now with the resonance actually placeable, which it\n' ...
         '  was not in cl_octave or twomap2.\n']);
disp('VENT_RES_DONE');
