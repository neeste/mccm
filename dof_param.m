function dof_param
% STAGE 3 CHECK: is dof now a parameter, and what does nch=1,dof=3 actually need?
%
% Two questions, deliberately separate:
%   (a) does resolve_dof honour an explicit request without breaking defaults?
%   (b) does the capstone workhorse configuration -- 1 chamber, 3 DOFs -- RUN?
% (b) is the plan's Stage 3 deliverable. If it does not run, what it needs is
% the useful output, not a pass/fail.

fprintf('\n  (a) resolve_dof: defaults unchanged, explicit request honoured\n');
fprintf('      config                        | dof | expected\n');
fprintf('      %s\n', repmat('-',1,52));
cases = { 'modpar26(1) default',        modpar26(1),                 2
          'modpar26(3) default (d3int)', modpar26(3),                3
          'modpar26(4) default (clvoct)', modpar26(4),               4 };
for i = 1:size(cases,1)
    pa = cases{i,2};
    if (isfield(pa,'dof')), pa = rmfield(pa,'dof'); end
    d = probe_dof(pa);
    fprintf('      %-29s | %3d | %d %s\n', cases{i,1}, d, cases{i,3}, tf(d==cases{i,3}));
end
p = modpar26(1); if (isfield(p,'dof')), p = rmfield(p,'dof'); end
p.dof = 3; d = probe_dof(p);
fprintf('      %-29s | %3d | %d %s   <== the capstone request\n', 'modpar26(1) + pa.dof=3', d, 3, tf(d==3));
p = modpar26(4); if (isfield(p,'dof')), p = rmfield(p,'dof'); end
p.dof = 2; d = probe_dof(p);
fprintf('      %-29s | %3d | %d %s   (floor protects m>=4)\n', 'modpar26(4) + pa.dof=2', d, 4, tf(d==4));

fprintf('\n  (b) does nch=1, dof=3 run?\n');
p = modpar26(1); if (isfield(p,'dof')), p = rmfield(p,'dof'); end
p.dof = 3; p.d3int = 1;
have5 = isfield(p,'k5o');
fprintf('      k5o/r5o/m5o present in the m=1 parameter set: %s\n', tf2(have5));
if (~have5)
    q = modpar26(3);
    p.k5o=q.k5o; p.k5e=q.k5e; p.k5q=q.k5q;
    p.r5o=q.r5o; p.r5e=q.r5e; p.r5q=q.r5q;
    p.m5o=q.m5o; p.m5e=q.m5e; p.m5q=q.m5q;
    fprintf('      -> borrowed from modpar26(3) for this probe\n');
end
p.isv = [1136 1005 840 655 466 273 80];
try
    evalc('S = tdm26(0,p,0,0);');
    fin = all(isfinite(S.d1(:)));
    fprintf('      RUNS: %s   max|d1| %.3e   d3 saved: %s\n', tf2(true), max(abs(S.d1(:))), ...
            tf2(isfield(S,'d3') && any(S.d3(:)~=0)));
    if (~fin), fprintf('      but NON-FINITE d1 -- unstable, not usable\n'); end
catch e
    fprintf('      DOES NOT RUN: %s\n', e.message);
end
fprintf(['\n  NOTE the d3 reference convention at m=1 is an OPEN question, not a\n' ...
         '  bug: macro_couple sets dref(3) only for m>=3, so at m=1 the third DOF\n' ...
         '  would be RELATIVE while the d3int force law (written for m=3b) assumes\n' ...
         '  ABSOLUTE. Flagged rather than guessed -- it is design-note decision 1\n' ...
         '  territory and changing it silently is how a refactor becomes a physics\n' ...
         '  change.\n']);
end

function d = probe_dof(pa)
% Via the CLICK path: dof_diag showed coupeig returns {lam,nch,boost,done,maxRe,
% maxRe_osc} with NO .pa, so the original probe caught an error for every case
% and reported -1 uniformly -- a probe failure that looked like a resolve_dof
% failure. The click path does return pa (vent_res.m relies on S.cp likewise).
pa.isv = [1136 840 466 80];
try, evalc('C = tdm26(0,pa,0,0);'); d = C.pa.dof; catch, d = -1; end
end
function s = tf(b),  if b, s='ok'; else, s='MISMATCH'; end, end
function s = tf2(b), if b, s='yes'; else, s='no'; end, end
