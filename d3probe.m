% D3PROBE -- is the third DOF LIVE at nch=3, and is d3int=1 fittable at all?
% Mistake guard #2/#3: identical output across a parameter change is plumbing,
% not physics. If k5/r5/m5 do not move the objective, fitting d3int=1 fits a
% 2-DOF model wearing a 3-DOF label, which is what every m=3b result on record
% turned out to be.
PROJ = '/Users/neely/Library/CloudStorage/OneDrive-FatherFlanagan''sBoysHome/STN/mod/mccm';
addpath(PROJ); cd(PROJ);
L=load('fit_nch3_surface.mat'); pa=L.R.pa;
if (isfield(pa,'hbmode')), pa=rmfield(pa,'hbmode'); end
p3=pa; p3.d3int=1;
fprintf('\n  k5o %g  r5o %g  m5o %g\n', gf(p3,'k5o'), gf(p3,'r5o'), gf(p3,'m5o'));
fprintf('\n=== does d3int=1 run, and is d3 LIVE? (fdm26 maperr) ===\n');
R2=fdm26(struct('pa',pa)); fprintf('  d3int=0            maperr %10.4f\n', R2.maperr);
R3=fdm26(struct('pa',p3)); fprintf('  d3int=1            maperr %10.4f\n', R3.maperr);
for f={'k5o','r5o','m5o'}
    q=p3; nm=f{1};
    if (~isfield(q,nm)), fprintf('  %-6s ABSENT from pa\n', nm); continue; end
    q.(nm)=q.(nm)*2;
    try
        Rq=fdm26(struct('pa',q));
        fprintf('  %-6s x2            maperr %10.4f   delta %+.3e %s\n', nm, Rq.maperr, ...
                Rq.maperr-R3.maperr, tern(abs(Rq.maperr-R3.maperr)<1e-9,'<-- INERT',''));
    catch e, fprintf('  %-6s x2 FAILED: %s\n', nm, e.message(1:min(40,end))); end
end
fprintf('\n=== time domain: does d3int=1 stay finite? ===\n');
try
    S=tdm26('wnr1',struct('fr',2,'lv',60,'pa',p3),0,0);
    fprintf('  wnr1 ok: tpk %.2f ms | d1mx %.3e | d3mx %.3e | d3 finite %d\n', ...
            S.tpk, S.dgn.d1mx, S.dgn.d3mx, S.dgn.d3fin);
catch e, fprintf('  wnr1 FAILED: %s\n', e.message(1:min(70,end))); end
fprintf('\n=== stability / amp at d3int=1 ===\n');
try
    evalc('C=tdm26(''coupeig'',struct(''pa'',p3));'); s=score26(p3,'fast',false);
    fprintf('  maxRe %+.2f | maxRe_osc %+.2f | amp %+.2f\n', C.maxRe, C.maxRe_osc, s.amp_gain);
catch e, fprintf('  FAILED: %s\n', e.message(1:min(70,end))); end
disp('D3PROBE_DONE');
function v=gf(s,f), if (isfield(s,f)), v=s.(f); else, v=NaN; end, end
function s=tern(c,a,b), if c, s=a; else, s=b; end, end
