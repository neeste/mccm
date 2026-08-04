% PVTEST -- the parameter mapping must not have moved under any saved fit.
PROJ = '/Users/neely/Library/CloudStorage/OneDrive-FatherFlanagan''sBoysHome/STN/mod/mccm';
addpath(PROJ); cd(PROJ);
H=parfit26('handles');
fits={'fit_nch3_surface.mat','fit_nch3_2dof.mat','fit_nch1_hbsc.mat','refit_c1broad.mat','refit_m4_full.mat'};
fprintf('\n  %-24s %6s %10s %14s\n','fit','nch','pv len','maperr delta');
for i=1:numel(fits)
    if (~exist(fits{i},'file')), continue; end
    L=load(fits{i}); if (isfield(L,'R')), p=L.R.pa; else, f=fieldnames(L); p=L.(f{1}); if isfield(p,'pa'), p=p.pa; end; end
    base=modpar26(p.m); if (isfield(p,'d3int')), base.d3int=p.d3int; end
    pv=H.getpar(p); pe=H.setpar(base,pv,[]);
    try
        a=fdm26(struct('pa',p)); b=fdm26(struct('pa',pe));
        fprintf('  %-24s %6d %10d %14.3e %s\n', fits{i}, p.m, numel(pv), abs(a.maperr-b.maperr), ...
                tern(abs(a.maperr-b.maperr)<1e-9,'','  <-- MAPPING MOVED'));
    catch e, fprintf('  %-24s ERROR %s\n', fits{i}, e.message(1:min(40,end))); end
end
fprintf('\n=== old (short) pv still loads and means the same? ===\n');
L=load('fit_nch3_surface.mat'); p=L.R.pa; nc=numel(p.chsz);
pv=H.getpar(p); pvold=pv(1:30+nc+2);          % simulate a pre-2026-08-04 vector
b1=H.setpar(modpar26(3),pv,[]); b2=H.setpar(modpar26(3),pvold,[]);
r1=fdm26(struct('pa',b1)); r2=fdm26(struct('pa',b2));
fprintf('  full pv maperr %.6f | truncated pv maperr %.6f | delta %.3e\n', ...
        r1.maperr, r2.maperr, abs(r1.maperr-r2.maperr));
fprintf('  (delta should be ~0: the appended entries equal modpar26(3) defaults here)\n');
fprintf('\n=== nch=1 must NOT gain a fabricated third DOF ===\n');
p1=modpar26(1); pv1=H.getpar(p1); q1=H.setpar(modpar26(1),pv1,[]);
fprintf('  k5o present after setpar at nch=1: %d (want 0)\n', isfield(q1,'k5o'));
disp('PVTEST_DONE');
function s=tern(c,a,b), if c, s=a; else, s=b; end, end
