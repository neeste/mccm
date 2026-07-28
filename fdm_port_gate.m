% GATE THE fdm26 PORT BEFORE TRUSTING ANY maperr IT PRODUCES.
%
% This project has FOUR previous failed fdm26 validation harnesses on record
% (reduced over frequency instead of place; guessed Ybm when the field is Yb;
% called a local subfunction externally; an fdsolve that ignored f). So the port
% gets a regression gate first and interpretation second.
%
% WHAT WAS PORTED
%   1. zv, a ninth imped output: the CL vent impedance kv/s + rv + (m5/clvent)*s.
%      The inertance m5/clvent is forced, not chosen -- tdm26's a2 stamp
%      G = clvent*m1/m5 is exactly the coefficient of a DOF of that mass.
%   2. nested: D(1,:) = [1 0 0 -1], moving d1 from ST<->SV to ST<->CL. Only row
%      1 changes. This must match tdm26's fold_p pickup (P_ST - P_CL).
%   3. clvent: a four-entry conductance stamp of Yv=1/zv on Y between CL and
%      pa.clvtgt -- the frequency-domain counterpart of tdm26's a2 stamp.
%
% ROW 1 IS THE REGRESSION GATE. Appended, no vent, stock areas must return
% maperr 1015.2 EXACTLY, the value clsweep2 reported in all ten of its rows
% before the port. Any drift there means the port changed the m=4 path it was
% supposed to leave alone, and every other row is void.
% ROW 2 is the same check at the carved areas, where the pre-port value was
% 1020.3 (carve_sv, vent_sv and vent_ss all agree on this).
%
% ROWS 3-6 are new capability. Note maperr was PINNED at the appended value in
% every nested row all session because fdm26 ignored both flags, so those
% numbers described a different model than tdm26 was amplifying. These are the
% first that describe the same model.
%
% bf_range comes from tdm26 and maperr from fdm26, so reading them side by side
% is itself a consistency check: if the port is wrong they will disagree about
% whether a configuration has a sane map.

A_STK = [0.95 0.05 1.00 0.05];
A_CRV = [0.95 0.05 0.95 0.05];

%        label                      chsz    nested vent  vtgt clvk
cfg = { 'REGRESSION appended stock', A_STK, 0, 0,    2, 0
        'REGRESSION appended carve', A_CRV, 0, 0,    2, 0
        'appended + vent 3 -> SS  ', A_CRV, 0, 3,    2, 0
        'nested sealed            ', A_CRV, 1, 0,    2, 0
        'nested + vent 3 -> SS    ', A_CRV, 1, 3,    2, 0
        'nested + vent 3 -> SV    ', A_CRV, 1, 3,    3, 0 };

fprintf('\n  PRE-PORT maperr: appended stock 1015.2 | appended carve 1020.3\n');
fprintf('  Rows 1-2 MUST reproduce those exactly or the port is void.\n');
fprintf('  m=3b reference maperr 499.3. Lower is a better CF map.\n\n');
fprintf('  config                    | maperr  | amp d1  | maxRe    | range mono\n');
fprintf('%s\n', repmat('-',1,74));
for i = 1:size(cfg,1)
    pa = modpar26(4); pa.chsz = cfg{i,2};
    pa.nested = cfg{i,3}; pa.clvent = cfg{i,4}; pa.clvtgt = cfg{i,5};
    if (cfg{i,6} ~= 0), pa.clvk = cfg{i,6}; end
    try
        S = score26(pa, 'fast', false);
    catch e
        fprintf('  %-25s | FAILED: %s\n', cfg{i,1}, e.message); continue
    end
    mk = '';
    if (i==1 && abs(S.maperr-1015.2) < 0.15), mk = ' <== GATE OK'; end
    if (i==2 && abs(S.maperr-1020.3) < 0.15), mk = ' <== GATE OK'; end
    if (i<=2 && isempty(mk)), mk = ' <== GATE FAILED'; end
    fprintf('  %-25s | %7.1f | %+6.2f | %+8.1f | %5.2f %-4s%s\n', ...
        cfg{i,1}, S.maperr, S.amp_gain, S.maxRe, S.bf_range, S.bf_mono, mk);
end
fprintf(['\n  If rows 1-2 hold, row 5 is the first honest map measurement of the\n' ...
         '  construction that met the amplifier bar, and the clean-map half of\n' ...
         '  SN''s bar becomes testable for the first time this session.\n']);
disp('FDM_PORT_GATE_DONE');
