% Verify the consolidation: modpar26(3) must now carry d3, and setting
% d3int=0 must reproduce the legacy 2-DOF model exactly.
p3 = modpar26(3);
fprintf('\nmodpar26(3): d3int=%d  dof-relevant fields present: k5o=%d m5o=%d\n', ...
        p3.d3int, isfield(p3,'k5o'), isfield(p3,'m5o'));
pl = p3; pl.d3int = 0;                       % legacy 2-DOF
ISV=[1136 1005 840 655 466 273 80];
p3.isv=ISV; pl.isv=ISV;
evalc('Sa=tdm26(0,p3,0,0);');   % consolidated (3 DOF)
evalc('Sb=tdm26(0,pl,0,0);');   % legacy (2 DOF)
fprintf('max |d1 consolidated - d1 legacy| = %.3e   (scale %.3e)\n', ...
        max(abs(Sa.d1(:)-Sb.d1(:))), max(abs(Sb.d1(:))));
fprintf('max |ped consolidated - ped legacy| = %.3e\n', max(abs(Sa.ped(:)-Sb.ped(:))));
fprintf('consolidated d3 (shadow RL) max = %.3e  ratio to d1 = %.2f\n', ...
        max(abs(Sa.d3(:))), max(abs(Sa.d3(:)))/max(abs(Sa.d1(:))));
fprintf('legacy d3 (must be all zero) max = %.3e\n', max(abs(Sb.d3(:))));
disp('CONSOL_VERIFY_DONE');
