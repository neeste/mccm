% Which response does the m=1 partition impedance actually drive?
% Code reading says fdm26.m:983 (zk(m,:)=1) destroys zk(1,1) before y is built
% at :1009, so y=1 and the BM response should be INSENSITIVE to k1/m1. The
% distillation pre-flight measured k1o/k1e as live. Both cannot be true.
% 50% perturbations -- unmissable if the parameter is live at all.
FL = 500*2.^(-1:5);
p  = modpar26(1);
R0 = fdm26(struct('cfmap',1,'pa',p,'flst',FL));
fprintf('\n  baseline xpk_bm: %s\n', num2str(R0.xpk_bm,'%.4f '));
fprintf('  baseline xpk_hb: %s\n\n', num2str(R0.xpk_hb,'%.4f '));
fprintf('  param      | max|d xpk_bm| | max|d xpk_hb|\n');
fprintf('  %s\n', repmat('-',1,48));
for f = {'k1o','m1o','k2o','k3o','k5o','m5o'}
    q = p; if (~isfield(q,f{1})), fprintf('  %-10s | (absent)\n', f{1}); continue; end
    q.(f{1}) = q.(f{1}) * 1.5;
    R = fdm26(struct('cfmap',1,'pa',q,'flst',FL));
    db = max(abs(R.xpk_bm - R0.xpk_bm));
    dh = max(abs(R.xpk_hb - R0.xpk_hb));
    fprintf('  %-10s | %13.3e | %13.3e\n', f{1}, db, dh);
end
fprintf(['\n  If xpk_bm is INSENSITIVE to k1o/m1o, the code reading is right and\n' ...
         '  m=1''s BM response does not carry its own partition impedance -- which\n' ...
         '  would make maperr 104.6, the best in the project, measure something\n' ...
         '  other than what we assume. That matters more than the distillation.\n']);
