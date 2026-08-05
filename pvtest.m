% PVTEST -- the parameter mapping must not have moved under any saved fit.
PROJ = '/Users/neely/Library/CloudStorage/OneDrive-FatherFlanagan''sBoysHome/STN/mod/mccm';
addpath(PROJ); cd(PROJ);
H=parfit26('handles');
fits={'fit_nch3_surface.mat','fit_nch3_2dof.mat','fit_nch1_hbsc.mat','refit_c1broad.mat','refit_m4_full.mat'};
% NAMES THE FIELD, not just the damage (2026-08-04). This test reported the m=4
% loss as a bare "136.5 <-- MAPPING MOVED" for as long as it has existed, which
% is a number nobody could act on: it says something is wrong without saying what,
% so the defect sat unfixed and got attributed to a guess (k5q/r5q/m5q) that
% turned out to be the wrong field entirely. It was clvent, the CL vent
% conductance, saved 0.5 against a default of 3. H.lost is the shipped check --
% NOT a paraphrase of it, which is this project's standing rule for tests.
% It is also a struct diff, so it catches fields that move the time domain
% without moving maperr, which the delta column is blind to by construction.
fprintf('\n  %-24s %6s %10s %14s  %s\n','fit','nch','pv len','maperr delta','fields lost');
for i=1:numel(fits)
    if (~exist(fits{i},'file')), continue; end
    L=load(fits{i}); if (isfield(L,'R')), p=L.R.pa; else, f=fieldnames(L); p=L.(f{1}); if isfield(p,'pa'), p=p.pa; end; end
    base=modpar26(p.m); if (isfield(p,'d3int')), base.d3int=p.d3int; end
    pv=H.getpar(p); pe=H.setpar(base,pv,[]);
    lost=H.lost(p,pe,[]);
    try
        a=fdm26(struct('pa',p)); b=fdm26(struct('pa',pe));
        fprintf('  %-24s %6d %10d %14.3e  %s%s\n', fits{i}, p.m, numel(pv), abs(a.maperr-b.maperr), ...
                tern(isempty(lost),'-',strjoin(lost,', ')), ...
                tern(abs(a.maperr-b.maperr)<1e-9,'','   <-- MAPPING MOVED'));
    catch e, fprintf('  %-24s ERROR %s\n', fits{i}, e.message(1:min(40,end))); end
end
fprintf('  (a named field here is FIXED by passing it through opts.pin; parfit26\n');
fprintf('   now REFUSES to start rather than resetting it silently.)\n');
fprintf('\n=== old (short) pv still loads and means the same? ===\n');
L=load('fit_nch3_surface.mat'); p=L.R.pa; nc=numel(p.chsz);
pv=H.getpar(p); pvold=pv(1:30+nc+2);          % simulate a pre-2026-08-04 vector
b1=H.setpar(modpar26(3),pv,[]); b2=H.setpar(modpar26(3),pvold,[]);
r1=fdm26(struct('pa',b1)); r2=fdm26(struct('pa',b2));
fprintf('  full pv maperr %.6f | truncated pv maperr %.6f | delta %.3e\n', ...
        r1.maperr, r2.maperr, abs(r1.maperr-r2.maperr));
fprintf('  (delta should be ~0: the appended entries equal modpar26(3) defaults here)\n');
fprintf('\n=== setpar must NOT FABRICATE a third DOF where the model has none ===\n');
% CORRECTED 2026-08-04. This asserted "k5o absent after setpar at nch=1", printed
% "1 (want 0)", and was a bug in the TEST, not the code: modpar26(1) SHIPS
% k5o/r5o/m5o (modpar26.m:9-11, seeded so "1 chamber, 3 DOFs" is reachable), so
% the field is present before setpar is ever called and its presence afterwards
% proves nothing. A test that prints a failure for correct code is worse than no
% test -- it trains the reader to ignore the line.
%
% The real property is that setpar_l CREATES nothing: it writes k5/r5/m5 only
% where the base already carries them (the isfield guards). So strip them from
% the base and check they stay gone, which is the case that would actually
% fabricate a DOF.
p1=modpar26(1); pv1=H.getpar(p1);
b0=rmfield(modpar26(1),{'k5o','r5o','m5o','k5e','r5e','m5e'});
q1=H.setpar(b0,pv1,[]);
made = intersect(fieldnames(q1), {'k5o','r5o','m5o','k5e','r5e','m5e'});
fprintf('  base stripped of k5/r5/m5; fields setpar re-created: %d (want 0) %s\n', ...
        numel(made), tern(isempty(made),'','  <-- FABRICATED A DOF'));
% and the ordinary path must still round-trip them where they DO exist
q2=H.setpar(modpar26(1),pv1,[]);
fprintf('  with the stock base, k5o round-trips: %.6g -> %.6g\n', p1.k5o, q2.k5o);
disp('PVTEST_DONE');
function s=tern(c,a,b), if c, s=a; else, s=b; end, end
