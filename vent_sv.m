% VENT CL TO SV, its parent, and see whether the amplifier comes back.
%
% WHERE THIS COMES FROM. The 2x2 showed the nested topology kills the amplifier
% outright: amp d1 -0.01 dB, maxRe exactly +0.0, and the sign fix worth +77 dB
% under the appended chain is worth 0.34 dB under nested. A force whose sign no
% longer matters is not reaching the travelling wave.
%
% THE MECHANISM. Nesting puts CL (area 0.05) in series in the d1 path, against
% ST at 0.95. The classical wave propagates by ST<->SV area exchange, and the
% nested chain replaces that direct path with three partitions in series whose
% middle link is 20x too small. The series inertance chokes the wave.
%
% WHY THE OLD VENT MADE IT WORSE RATHER THAN BETTER. It stamped G onto the SAME
% (ST,CL) matrix entries that already carry d1, so it stiffened the BM's own
% pressure difference instead of opening a second path. It attacked the drive.
%
% THE FIX. CL is carved from SV, so it should communicate with SV. Venting to
% the parent restores an ST<->SV exchange route THROUGH CL while leaving
% p_ST - p_CL intact as the BM drive. Large clvent should then recover
% appended-like amplification, which is this knob's reduction gate. It will not
% land exactly on appended even at infinite vent, since SS still sits between CL
% and SV in the chain.
%
% GATE ROWS. vent=0 must reproduce the sealed nested row from the 2x2 exactly
% (-0.01 / +0.0), confirming the refactor to an additive vent changed nothing.
% The vent->ST row must still fail, confirming the direction is what matters and
% not merely the presence of a vent.
%
% NOTE clcouple is inert under nested: that branch never reads cc. So the m=4
% reduction gate does not currently exist in the nested form.

A_CRV = [0.95 0.05 0.95 0.05];   % CL carved from SV, total conserved at 2.00
ST = 1; SV = 3;                  % vent targets, as chamber indices

cfg = { 'appended (reference)', 0, 0,    SV
        'nested sealed  [gate]', 1, 0,    SV
        'nested vent 0.1 -> SV', 1, 0.1,  SV
        'nested vent 0.5 -> SV', 1, 0.5,  SV
        'nested vent 1   -> SV', 1, 1,    SV
        'nested vent 3   -> SV', 1, 3,    SV
        'nested vent 10  -> SV', 1, 10,   SV
        'nested vent 30  -> SV', 1, 30,   SV
        'nested vent 100 -> SV', 1, 100,  SV
        'nested vent 10  -> ST', 1, 10,   ST };

fprintf('\n  TARGETS  m=3b: amp d1 +81.15 | maxRe +19.3 | maperr 499.3\n');
fprintf('  2x2 rows appended SV=0.95: +80.62 | maxRe +4524.8 | nested sealed: -0.01 | +0.0\n');
fprintf('  SN bar for m=4: clean map AND amp 40-60 dB\n\n');
fprintf('  config                | amp d1  amp d2  amp d3 | maxRe     | range mono | maperr\n');
fprintf('%s\n', repmat('-',1,88));
R = struct('nm',{},'S',{});
for i = 1:size(cfg,1)
    pa = modpar26(4); pa.chsz = A_CRV;
    pa.nested = cfg{i,2}; pa.clvent = cfg{i,3}; pa.clvtgt = cfg{i,4};
    try
        S = score26(pa, 'fast', false);
    catch e
        fprintf('  %-21s | FAILED: %s\n', cfg{i,1}, e.message); continue
    end
    flag = '';
    if (S.maxRe > 24), flag = ' UNSTABLE'; end
    if (S.amp_gain>=40 && S.amp_gain<=60 && strcmp(S.bf_mono,'ok') && S.maxRe<=24)
        flag = ' <== MEETS THE BAR';
    end
    fprintf('  %-21s | %+6.2f %+7.2f %+7.2f | %+9.1f | %5.2f %-4s | %6.1f%s\n', ...
        cfg{i,1}, S.amp_gain, S.amp_d2, S.amp_d3, S.maxRe, S.bf_range, ...
        S.bf_mono, S.maperr, flag);
    R(end+1).nm = cfg{i,1}; R(end).S = S; %#ok<SAGROW>
    save('vent_sv.mat','R');
end
fprintf(['\n  READ: does amp d1 climb off zero as the vent to SV opens? If it\n' ...
         '  recovers toward the appended +80 while maxRe stays in the healthy\n' ...
         '  band, the nested topology is viable and the vent direction was the\n' ...
         '  whole fault. If it stays flat at zero for every vent strength, then\n' ...
         '  routing the wave through a 0.05-area chamber cannot be rescued by a\n' ...
         '  shunt, and CL should be a side space off the partition rather than a\n' ...
         '  series element in the d1 path.\n']);
disp('VENT_SV_DONE');
