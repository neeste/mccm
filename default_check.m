% VERIFY THE NEW modpar26c4 DEFAULT.
%
% modpar26c4 now ships the nested construction: chsz carved from SV (total
% conserved at 2.00), nested=1, clvent=3.0 to SS, clvoct=4.0.
%
% clvoct is a NEW parameter resolved inside tdm26 and fdm26 rather than in the
% parameter file, because the stiffness it implies is place-dependent (it tracks
% k1/m1/m5); a scalar clvk in modpar26c4 would be correct at exactly one place.
% Both solvers must resolve it IDENTICALLY or maperr describes a different model
% than the amplifier measurement, which is the exact fault the fdm26 port was
% built to remove.
%
% ROW 1 is the check that matters: the bare default must reproduce the measured
% clvoct-4 row (maperr 525.2, amp +56.59, maxRe +2.3, contrast 6.4). Those came
% from an explicit place-dependent clvk computed in vent_res2, so agreement
% proves the in-solver clvoct resolution matches the hand-computed stiffness.
%
% ROW 2 checks the dof promotion. If clvoct did not trigger it, the vent would
% carry a resonance the time march has no state for and row 1 would silently
% return the clvk=0 numbers (522.7 / +55.08 / 4.2) -- a failure that reads as
% success. Row 2 IS that comparison, made explicit.
%
% ROW 3 confirms the documented escape hatch restores the legacy model.
%
% NOTE built without local functions or anonymous handles. Local functions in a
% SCRIPT get their own workspace, and anonymous functions calling them is a
% further hazard; this cost a run earlier today in cl_octave2.

fprintf('\n  TARGET default (= measured clvoct 4.0 row):\n');
fprintf('    maperr 525.2 | amp +56.59 | maxRe +2.3 | contrast 6.4\n');
fprintf('  If clvoct silently failed, expect the clvk=0 row instead:\n');
fprintf('    maperr 522.7 | amp +55.08 | maxRe +3.9 | contrast 4.2\n');
fprintf('  Legacy appended: maperr 1020.3 | amp +80.62 | maxRe +4524.8\n\n');

p1 = modpar26(4);                                  % the new default
fprintf('  chsz      = [%s]  (sum %.2f, must be 2.00)\n', ...
    strtrim(sprintf('%.2f ', p1.chsz)), sum(p1.chsz));
fprintf('  nested %d | clvtgt %d | clvent %.1f | clvoct %.1f\n', ...
    p1.nested, p1.clvtgt, p1.clvent, p1.clvoct);
S0 = tdm26(0, p1, 0, 0);
fprintf('  pa.dof resolved to %d (MUST be 4, else clvoct has no state)\n\n', S0.pa.dof);

p2 = modpar26(4);                                  % same, resonance disabled
p2 = rmfield(p2, 'clvoct');

p3 = modpar26(4);                                  % documented legacy recovery
p3.chsz = [0.95 0.05 1.0 0.05];
p3.nested = 0;
for fn = {'clvent','clvoct','clvtgt'}
    if (isfield(p3, fn{1})), p3 = rmfield(p3, fn{1}); end
end

lbl = {'new default        ', 'default, clvoct off', 'legacy appended    '};
PA  = {p1, p2, p3};
fprintf('  config              | maperr  | amp d1  | maxRe    | contrast | range mono\n');
fprintf('%s\n', repmat('-',1,80));
for i = 1:numel(PA)
    try
        S = score26(PA{i}, 'fast', false);
    catch e
        fprintf('  %-19s | FAILED: %s\n', lbl{i}, e.message); continue
    end
    fprintf('  %-19s | %7.1f | %+6.2f | %+8.1f | %8.1f | %5.2f %-4s\n', ...
        lbl{i}, S.maperr, S.amp_gain, S.maxRe, S.contrast, S.bf_range, S.bf_mono);
end
fprintf(['\n  Row 1 must hit 525.2 / +56.59 / 6.4. Landing on row 2''s numbers\n' ...
         '  instead means clvoct resolved to nothing and the default is not what\n' ...
         '  it claims. Row 3 must return the legacy appended numbers, confirming\n' ...
         '  the escape hatch documented in modpar26c4 actually works.\n']);
disp('DEFAULT_CHECK_DONE');
