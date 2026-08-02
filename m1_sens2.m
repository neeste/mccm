% Does the third DOF activate at dof=3, and stay OUT of the way at dof=2?
% m1_sens.m tested only modpar26(1) with no pa.dof, so use3 was false and the
% 2-DOF path ran -- correct behaviour, wrong test. Both branches now.
FL = 500*2.^(-1:5);
for d = [2 3]
    p = modpar26(1); p.dof = d;
    R0 = fdm26(struct('cfmap',1,'pa',p,'flst',FL));
    fprintf('\n  === pa.dof = %d ===\n', d);
    fprintf('  xpk_bm: %s\n', num2str(R0.xpk_bm,'%.4f '));
    fprintf('  param      | max|d xpk_bm| | max|d xpk_hb| | live\n');
    fprintf('  %s\n', repmat('-',1,56));
    for f = {'k1o','k5o','r5o','m5o'}
        q = p; if (~isfield(q,f{1})), continue; end
        q.(f{1}) = q.(f{1}) * 1.5;
        R = fdm26(struct('cfmap',1,'pa',q,'flst',FL));
        db = max(abs(R.xpk_bm-R0.xpk_bm)); dh = max(abs(R.xpk_hb-R0.xpk_hb));
        lv = 'no '; if (max(db,dh) > 1e-12), lv = 'YES'; end
        fprintf('  %-10s | %13.3e | %13.3e | %s\n', f{1}, db, dh, lv);
    end
end
% dof=2 must reproduce the pre-edit numbers exactly: the third DOF is opt-in.
p2 = modpar26(1); p2.dof = 2; R2 = fdm26(struct('cfmap',1,'pa',p2,'flst',FL));
p0 = modpar26(1);                R0 = fdm26(struct('cfmap',1,'pa',p0,'flst',FL));
fprintf('\n  dof=2 vs no-dof-field: max|d xpk_bm| = %.3e  (must be 0)\n', ...
        max(abs(R2.xpk_bm-R0.xpk_bm)));
fprintf(['\n  READ: k5o/r5o/m5o must be INERT at dof=2 (opt-in preserved) and LIVE\n' ...
         '  at dof=3 (the blocker removed). Anything else means the gate is wrong.\n']);
