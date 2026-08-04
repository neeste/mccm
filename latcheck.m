% LATCHECK -- acceptance tests for the soft-argmax latency detector.
%
% 1. LEGACY EXACT. latsoft=Inf must reproduce the numbers drvcheck recorded
%    under the old code (slope 0.4118 / 0.5058 / lvl_c 6.399). If not, the edit
%    changed something it should not have.
% 2. CLEAN TRACES UNPERTURBED. stock nch=1 has shoulder 0.0000 over all 16
%    cells, so the soft-argmax must return the SAME latencies as the argmax.
%    This is the test that the smoothing does not tax a healthy model.
% 3. DRIVE DEPENDENCE SHRINKS. The whole point: |slope(bm)-slope(default)| was
%    0.0940, i.e. 78x the gap c3fit spent 17 h closing.
% 4. CONTINUITY. The hbsc axis was non-monotonic under the old detector
%    (0.4151 -> 0.3665 -> 0.4379 -> 0.4452 across 0.0406/0.06/0.10/0.20) and I
%    wrongly read that as noise. Under a continuous detector it should be
%    smooth. THIS IS THE REAL ACCEPTANCE TEST -- the others only check I did
%    not break anything.
PROJ = '/Users/neely/Library/CloudStorage/OneDrive-FatherFlanagan''sBoysHome/STN/mod/mccm';
addpath(PROJ);
S3=load('/Users/neely/mccm_runs/c3fit_seg3.mat'); pa=S3.R.pa;
nod=@(p) rmf(p,'hbmode');

fprintf('\n=== 1. LEGACY EXACT (latsoft=Inf vs drvcheck under the OLD code) ===\n');
pL=nod(pa); pL.latsoft=Inf;              mL=abr_metric(pL,false);
pLb=pa; pLb.hbmode='bm'; pLb.latsoft=Inf; mLb=abr_metric(pLb,false);
fprintf('  default  slope %.4f (drvcheck 0.4118)  lvl_c %.3f (6.399)\n', mL.slope, mL.level_c);
fprintf('  bm       slope %.4f (drvcheck 0.5058)  lvl_c %.3f (6.244)\n', mLb.slope, mLb.level_c);
okleg = abs(mL.slope-0.4118)<5e-4 && abs(mLb.slope-0.5058)<5e-4;
fprintf('  %s\n', ternx(okleg,'LEGACY PATH EXACT','*** LEGACY PATH CHANGED -- investigate ***'));

fprintf('\n=== 2. CLEAN TRACES UNPERTURBED (stock nch=1, shoulder 0) ===\n');
p1=modpar26(1); p1.latsoft=Inf;  m1L=abr_metric(p1,false);
p1s=modpar26(1); p1s.latsoft=8;  m1S=abr_metric(p1s,false);
dlat=max(abs(m1L.lat(:)-m1S.lat(:)));
fprintf('  shoulder %.4f | max |lat_soft - lat_argmax| = %.4f ms | slope %.4f vs %.4f\n', ...
        m1L.shoulder, dlat, m1S.slope, m1L.slope);
fprintf('  %s\n', ternx(dlat<1e-9,'IDENTICAL on a clean single-peaked WNR', ...
                        'soft-argmax perturbs a clean trace -- softp may be too low'));

fprintf('\n=== 3. DRIVE DEPENDENCE (softp=8) ===\n');
pd=nod(pa); pd.latsoft=8;                mD=abr_metric(pd,false);
pb=pa; pb.hbmode='bm'; pb.latsoft=8;     mB=abr_metric(pb,false);
old=0.0940; new=abs(mB.slope-mD.slope);
fprintf('  default  slope %.4f  lvl_c %.3f  shoulder %.4f\n', mD.slope, mD.level_c, mD.shoulder);
fprintf('  bm       slope %.4f  lvl_c %.3f  shoulder %.4f\n', mB.slope, mB.level_c, mB.shoulder);
fprintf('  |bm - default| = %.4f   (was %.4f, %.1fx reduction)\n', new, old, old/max(new,1e-9));

fprintf('\n=== 4. CONTINUITY of the hbsc axis ===\n');
hb=[0.0406 0.06 0.10 0.20];
oldslope=[0.4151 0.3665 0.4379 0.4452];   % lcprobe, old detector
fprintf('  %8s %12s %12s\n','hbsc','slope OLD','slope NEW');
ns=zeros(size(hb));
for i=1:numel(hb)
    p=nod(pa); p.latsoft=8; p.hbsc=hb(i);
    m=abr_metric(p,false); ns(i)=m.slope;
    fprintf('  %8.4f %12.4f %12.4f\n', hb(i), oldslope(i), ns(i));
end
dold=diff(oldslope); dnew=diff(ns);
fprintf('  sign changes in successive differences:  OLD %d   NEW %d\n', ...
        sum(diff(sign(dold))~=0), sum(diff(sign(dnew))~=0));
fprintf('  %s\n', ternx(sum(diff(sign(dnew))~=0)==0, ...
        'MONOTONIC under the new detector -- the jaggedness was the peak switch', ...
        'still non-monotonic -- softp may be too high, or there is a second cause'));
save('/Users/neely/mccm_runs/latcheck.mat','mL','mLb','mD','mB','m1L','m1S','ns','hb');
disp('LATCHECK_DONE');

function p=rmf(p,f), if (isfield(p,f)), p=rmfield(p,f); end, end
function s=ternx(c,a,b), if c, s=a; else, s=b; end, end
