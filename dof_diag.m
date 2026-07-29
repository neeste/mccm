function dof_diag
% Diagnostic: (1) what does coupeig actually return? (2) where exactly does
% nch=1,dof=3 fail? No try/catch on the second -- the stack is the point.
fprintf('\n  (1) coupeig return fields (my probe assumed .pa and was wrong)\n');
S = tdm26('coupeig', struct('pa', modpar26(3)));
fprintf('      %s\n', strjoin(fieldnames(S)', ', '));
fprintf('      has .pa: %d\n', isfield(S,'pa'));

fprintf('\n  (2) does the CLICK path return pa? (vent_res.m used S.cp)\n');
p = modpar26(1); p.isv = [1136 840 466 80];
evalc('C = tdm26(0,p,0,0);');
fprintf('      fields: %s\n', strjoin(fieldnames(C)', ', '));
if (isfield(C,'pa')), fprintf('      C.pa.dof = %d  (m=1 default, expect 2)\n', C.pa.dof); end

fprintf('\n  (3) nch=1, dof=3 -- full stack, uncaught\n');
q = modpar26(3);
p = modpar26(1); p.dof = 3; p.d3int = 1;
p.k5o=q.k5o; p.k5e=q.k5e; p.k5q=q.k5q;
p.r5o=q.r5o; p.r5e=q.r5e; p.r5q=q.r5q;
p.m5o=q.m5o; p.m5e=q.m5e; p.m5q=q.m5q;
p.isv = [1136 840 466 80];
tdm26(0,p,0,0);
end
