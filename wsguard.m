% WSGUARD -- regression test for parfit26's WARM-START INTEGRITY guard.
%
% The guard refuses to start a fit whose warm start would be silently altered by
% setpar_l. The design risk is FALSE POSITIVES: a guard that also fired on
% legitimate warm starts would block real work, and would be turned off within a
% week. So the test that matters is not "does it catch the m=4 defect" -- it is
% "does it stay silent on everything that is fine". Both are asserted here.
%
% Zero model evaluations by construction: the check is a struct diff, which is
% why it can run before every fit at no cost.
PROJ = '/Users/neely/Library/CloudStorage/OneDrive-FatherFlanagan''sBoysHome/STN/mod/mccm';
addpath(PROJ); cd(PROJ);
fprintf('\n  testing: %s\n', which('parfit26'));
H = parfit26('handles');

fits = {'fit_nch3_surface.mat','fit_nch3_2dof.mat','fit_nch1_hbsc.mat', ...
        'refit_c1broad.mat','refit_m4_full.mat'};
fprintf('\n  %-24s %4s %-40s\n','fit','m','fields that would be silently reset');
for i=1:numel(fits)
    if (~exist(fits{i},'file')), continue; end
    L=load(fits{i});
    if (isfield(L,'R')), p=L.R.pa; else, f=fieldnames(L); p=L.(f{1}); if isfield(p,'pa'), p=p.pa; end; end
    base=modpar26(p.m); if (isfield(p,'d3int')), base.d3int=p.d3int; end
    rt = H.setpar(base, H.getpar(p), []);
    lost = H.lost(p, rt, []);
    if (isempty(lost)), s='(none -- survives)'; else, s=strjoin(lost,', '); end
    fprintf('  %-24s %4d %-40s\n', fits{i}, p.m, s);
end

% ---- the fix works: pinning clvent makes the m=4 fit survive ----
fprintf('\n  does opts.pin repair it?\n');
L=load('refit_m4_full.mat'); p=L.R.pa;
b1=modpar26(4);
fprintf('    unpinned      : %s\n', strjoin_e(H.lost(p, H.setpar(b1,H.getpar(p),[]), [])));
b2=modpar26(4); b2.clvent=p.clvent;                 % what opts.pin does
fprintf('    pin clvent    : %s\n', strjoin_e(H.lost(p, H.setpar(b2,H.getpar(p),[]), [])));

% ---- it must also catch the DOCUMENTED gampro trap ----
fprintf('\n  catches a warm start carrying gampro with no matching gpgrade?\n');
p3=modpar26(3); p3.gampro = exp(0.7*(((0:p3.n-1)')/(p3.n-1)-0.5));  % non-uniform, no gpgrade
lost3 = H.lost(p3, H.setpar(modpar26(3), H.getpar(p3), []), []);
fprintf('    %s %s\n', strjoin_e(lost3), tern(any(strcmp(lost3,'gampro')),'<== caught','<== MISSED'));

% ---- chszderive must NOT be flagged as a defect ----
fprintf('\n  chsz(3) derived is a constraint, not a loss?\n');
p4=modpar26(3); p4.chsz=[0.7 0.4 0.9];
lost4 = H.lost(p4, H.setpar(modpar26(3), H.getpar(p4), 3), 3);
fprintf('    %s %s\n', strjoin_e(lost4), tern(~any(strcmp(lost4,'chsz')),'<== correctly ignored','<== FALSE POSITIVE'));
disp('WSGUARD_DONE');
function s=strjoin_e(c), if isempty(c), s='(none)'; else, s=strjoin(c,', '); end, end
function s=tern(c,a,b), if c, s=a; else, s=b; end, end
