% SLOPEATLAS -- what does each fitted parameter cost in STABILITY per unit of
% ABR SLOPE it buys?
%
% WHY. The nch=1 objective is now entirely the slope term (map/gain/osc all
% exactly 0.0000 at the 90.97 optimum), and 150 evaluations left it at 0.606
% against a 0.413 target. The obvious lever, gampro, turned out to be CF-map-free
% but stability-bankrupt: ~2600 dB of oscillatory margin per unit of slope, so
% the entire 56.6 dB margin buys 0.021 of the 0.193 needed. Before spending a
% 12-20 h window on any fit, measure the SAME exchange rate for every parameter
% already in the vector. The fit cannot do better than its best-rate parameter.
%
% WHAT IT REPORTS, per parameter, from central differences at +-2%:
%   dslope  d(slope)/d(ln p)     how much slope a relative change buys
%   dosc    d(maxRe_osc)/d(ln p) what it costs in oscillatory margin
%   dmap    d(maperr)/d(ln p)    whether it also wrecks the CF map
%   RATE    |dslope| / |dosc|    slope per dB of margin -- THE COLUMN THAT MATTERS
%
% Read it as: with M dB of margin in hand, the reachable slope change is
% RATE*M. Baseline margin is 56.6 dB and the needed slope change is 0.193, so a
% parameter is only viable if RATE > 0.193/56.6 = 3.4e-3. gampro scores 3.8e-4,
% an order of magnitude short -- that threshold is the point of the exercise.
%
% Central differences, not one-sided: the gampro cliff (54.8 dB of margin for a
% 1% step) is exactly the shape that makes a one-sided estimate a coin flip.
t0=tic;
L=load('sweep_nch1b.mat'); pa0=L.R.pa;
H=parfit26('handles'); nm=H.parnames(); pv0=H.getpar(pa0);
nc=numel(pa0.chsz); gpi=30+nc+1;
idx=[1 2 3 4 5 7 8 9 10 11 13 20 21 31 32 gpi];   % nch=1 default fitidx + grade
lbl=cell(1,numel(idx));
for j=1:numel(idx)
    i=idx(j);
    if     (i<=30),      lbl{j}=nm{i};
    elseif (i<=30+nc),   lbl{j}=sprintf('chsz(%d)',i-30);
    else,                lbl{j}='gampro-g';
    end
end
rel=0.02;                                   % +-2%
ev=@(pv) evalpv(pv,H,pa0);   % local fn at EOF; a script cannot share a workspace
                             % with a nested function, so pass H/pa0 explicitly

fprintf('\n== slope-vs-stability atlas, nch=1, from sweep_nch1b ==\n');
[sl0,mp0,os0]=ev(pv0);
fprintf('   baseline: slope %.4f  maperr %.2f  maxRe_osc %.1f\n', sl0, mp0, os0);
fprintf('   need dslope = %.3f, have margin = %.1f dB  ->  VIABLE needs RATE > %.1e\n\n', ...
        0.413-sl0, -os0, abs(0.413-sl0)/abs(os0));
fprintf('   param       dslope/dlnp   dosc/dlnp   dmap/dlnp |     RATE    verdict\n');
fprintf('   ---------------------------------------------------------------------\n');
need=abs(0.413-sl0); marg=abs(os0); thr=need/marg;
A=struct('name',lbl,'dslope',NaN,'dosc',NaN,'dmap',NaN,'rate',NaN);
for j=1:numel(idx)
    i=idx(j); p=pv0(i);
    if (i==gpi), hstep=rel; else, hstep=rel*abs(p); end   % grade is 0: absolute step
    if (hstep==0), hstep=rel; end
    pvp=pv0; pvp(i)=p+hstep;  [slp,mpp,osp]=ev(pvp);
    pvm=pv0; pvm(i)=p-hstep;  [slm,mpm,osm]=ev(pvm);
    % per unit RELATIVE change, so parameters of wildly different units compare
    den=2*hstep/max(abs(p),hstep);
    ds=(slp-slm)/den; do=(osp-osm)/den; dm=(mpp-mpm)/den;
    rt=abs(ds)/max(abs(do),eps);
    A(j).dslope=ds; A(j).dosc=do; A(j).dmap=dm; A(j).rate=rt;
    v='--';
    if (isfinite(rt))
        if (rt>thr), v='VIABLE'; else, v=sprintf('%.0fx short',thr/max(rt,eps)); end
    end
    if (~isfinite(ds)), v='eval failed'; end
    fprintf('   %-10s  %+11.4f  %+10.1f  %+10.2f | %8.2e  %s\n', lbl{j}, ds, do, dm, rt, v);
end
save('slopeatlas.mat','A','idx','lbl','sl0','mp0','os0','thr');
fprintf('\n   %.1f min.  Best RATE = the only parameter a long fit can ride.\n', toc(t0)/60);
disp('SLOPEATLAS_DONE');

function [sl,mp,os]=evalpv(pv,H,pa0)
sl=NaN; mp=NaN; os=NaN;
pa=H.setpar(pa0,pv,[]);
try, m=abr_metric(pa,false); if (m.ok), sl=m.slope; end, catch, end
try, S=score26(pa,'fast',false); mp=S.maperr; os=S.maxRe_osc; catch, end
end
