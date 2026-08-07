% HBCHECK -- gate + measurement for pa.hbrl under rlsplit.
% GATE: hbrl off must be bit-identical (the drive is d2 either way).
% MEASUREMENT: does the RL-driven MET restore the shear the topology removed?
PROJ='/Users/neely/Library/CloudStorage/OneDrive-FatherFlanagan''sBoysHome/STN/mod/mccm';
addpath(PROJ); cd(PROJ);
L=load('/Users/neely/mccm_runs/rlfit_seg2.mat'); pa=L.R.pa;
if (isfield(pa,'hbmode')), pa=rmfield(pa,'hbmode'); end
pa.d3int=1; pa.rlsplit=1;
p0=pa; p0.hbrl=0;  p1=pa; p1.hbrl=1;
fprintf('\n  %-16s %10s %10s %10s %12s %10s\n','hbrl','d1mx','d2mx','dhb/d1','ohcBM','amp dB');
for i=1:2
    q = {p0,p1}; q=q{i};
    try
        S=tdm26('wnr1',struct('fr',2,'lv',60,'pa',q),0,0); d=S.dgn;
        s=score26(q,'fast',false); R=fdm26(struct('pa',q));
        fprintf('  %-16d %10.3e %10.3e %10.4f %12.3e %+10.2f  maperr %.2f\n', ...
                q.hbrl, d.d1mx, d.d2mx, d.d2mx/max(d.d1mx,eps), d.ohcBM, s.amp_gain, R.maperr);
    catch e, fprintf('  hbrl=%d FAILED: %s\n', q.hbrl, e.message(1:min(60,end))); end
end
% GATE: off must reproduce the committed rlfit result exactly
R0=fdm26(struct('pa',p0)); Rr=fdm26(struct('pa',pa));
fprintf('\n  GATE hbrl=0 vs as-fitted: maperr delta %.3e %s\n', abs(R0.maperr-Rr.maperr), ...
        tern(abs(R0.maperr-Rr.maperr)<1e-9,'BIT-IDENTICAL','  <-- NOT A NO-OP, investigate'));
disp('HBCHECK_DONE');
function s=tern(c,a,b), if c, s=a; else, s=b; end, end
