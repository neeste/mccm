% LATATLAS -- what actually sets the nch=3 forward latency?
%
% THE QUESTION, POSED PROPERLY. hbsc is not inert because it is weak: swept, it
% moves slope 0.6506 -> 0.3665 across hbsc 0.005-0.06, an enormous lever. It is
% inert because it drags level_c with it nearly one-for-one (1.802 -> 12.499
% over the same range), and every objective used so far constrains level_c --
% explicitly (wlcb band) in c3fit, implicitly (the level axis of the surface) in
% s3fit. Two long runs left it at 0.0403 for that reason, not for lack of power.
%
% So the useful quantity is NOT d(slope)/d(param). It is the SELECTIVITY
%
%       S = |d slope| / |d level_c|      per identical relative perturbation
%
% A lever with S ~ 1 is locked: it buys slope and pays for it in level
% dependence. A lever with S >> 1 moves latency while leaving the level axis
% alone, which is what the objective has never been able to find.
%
% Also tracked, because a lever that works by re-introducing the double peak is
% worthless (that is what the whole detector episode was about):
%       d shoulder   -- must not rise
%
% METHOD WARNING carried from abr-tuning-levers, and it cost time before: a
% +-x% central difference is NOT trustworthy alone here. It rated ace VIABLE for
% the right reason by a worthless argument. This is a SCREEN to rank candidates;
% every surviving row must then be confirmed by a direct sweep over the range
% the claim needs. Nothing here is a result on its own.
%
% Scored from fit_nch3_surface.mat, the honestly-scored current best (shear
% drive, soft-argmax detector, shoulder 0.2612).
PROJ = '/Users/neely/Library/CloudStorage/OneDrive-FatherFlanagan''sBoysHome/STN/mod/mccm';
addpath(PROJ); cd(PROJ);
L=load('fit_nch3_surface.mat'); pa=L.R.pa;
if (isfield(pa,'hbmode')), pa=rmfield(pa,'hbmode'); end
H=parfit26('handles'); nm=H.parnames();
pv0=H.getpar(pa); np=numel(pv0);
lbl=cell(1,np);
for i=1:30, lbl{i}=nm{i}; end
nc=numel(pa.chsz);
for i=1:nc, lbl{30+i}=sprintf('chsz%d',i); end
lbl{30+nc+1}='gpgrade'; lbl{30+nc+2}='hbsc';

REL=0.05;
m0=abr_metric(pa,false);
fprintf('\n  BASE  slope %.4f  level_c %.3f  shoulder %.4f\n', m0.slope, m0.level_c, m0.shoulder);
fprintf('  perturbation +-%.0f%%, central difference, %d parameters\n\n', 100*REL, np);

ds=nan(1,np); dl=nan(1,np); dh=nan(1,np);
for i=1:np
    v=pv0(i);
    if (abs(v)>1e-12), dv=abs(v)*REL; else, dv=REL; end
    s=nan(1,2); l=nan(1,2); h=nan(1,2); ok=true;
    for k=1:2
        pvk=pv0; pvk(i)=v+(2*k-3)*dv;
        pk=H.setpar(pa,pvk,[]);
        try
            mk=abr_metric(pk,false);
            if (~mk.ok), ok=false; break; end
            s(k)=mk.slope; l(k)=mk.level_c; h(k)=mk.shoulder;
        catch
            ok=false; break;
        end
    end
    if (~ok || any(~isfinite([s l h]))), fprintf('  %-8s  unusable\n', lbl{i}); continue; end
    ds(i)=(s(2)-s(1))/2; dl(i)=(l(2)-l(1))/2; dh(i)=(h(2)-h(1))/2;
    fprintf('  %-8s  dslope %+8.5f  dlvlc %+8.4f  dshldr %+7.4f\n', lbl{i}, ds(i), dl(i), dh(i));
end

sel = abs(ds)./max(abs(dl),1e-9);
fprintf('\n=== RANKED BY SELECTIVITY  S = |dslope| / |dlevel_c| ===\n');
fprintf('  (S ~ 1 = locked, buys slope and pays in level dependence)\n');
fprintf('  (need |dslope| meaningful too -- a tiny lever with high S is useless)\n\n');
fprintf('  %-8s %10s %10s %10s %10s\n','param','S','dslope','dlvlc','dshldr');
[~,ix]=sort(sel,'descend','MissingPlacement','last');
for k=1:min(14,np)
    i=ix(k);
    if (~isfinite(sel(i))), continue; end
    fprintf('  %-8s %10.2f %+10.5f %+10.4f %+10.4f\n', lbl{i}, sel(i), ds(i), dl(i), dh(i));
end
fprintf('\n  hbsc for reference: S = %.2f\n', sel(30+nc+2));
save('/Users/neely/mccm_runs/latatlas.mat','ds','dl','dh','sel','lbl','pv0','m0');
disp('LATATLAS_DONE');
