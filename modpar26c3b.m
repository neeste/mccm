function pa = modpar26c3b
%MODPAR26C3B  Alias for modpar26(3), kept for scripts that reference it.
%
% As of the 2026-07-26 consolidation, m=3 CARRIES the internal third DOF by
% default, so modpar26c3b and modpar26(3) are the same model. The separate
% "m=3b" rung existed while the staging step was being validated; its identity
% gate passed exactly (d1, d2, ped bit-identical to the 2-DOF m=3), which is
% what justified folding it in.
%
% Use pa.d3int = 0 for the legacy 2-DOF behaviour.
pa = modpar26(3);
end
