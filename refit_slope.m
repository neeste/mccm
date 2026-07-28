function R = refit_slope(fitidx, wslope, maxfe)
% REFIT_SLOPE  Month-2 JOINT re-fit: minimize tuning/MAP error AND the ABR
%   group-delay slope error together, over a subset of the 30 impedance params
%   (getpar/setpar order), then verify sub-criticality (coupeig).
%
%     J(pv) = maperr(pv) + wslope * |d(pv) - 0.413|
%
%   The single-param search showed slope and tuning fight along any one axis;
%   here the optimizer trades multiple coefficients at once to satisfy both.
%
%   fitidx : indices into the 30-param vector to fit (default: BM block
%            [k1o r1o m1o k1e r1e m1e] = [1 2 3 11 12 13]).
%   wslope : weight on the slope term (default 500; maperr~100, |d-.413|~0.28).
%   maxfe  : fminsearch MaxFunEvals (default 200).

if (nargin<1 || isempty(fitidx)), fitidx=[1 2 3 11 12 13]; end
if (nargin<2 || isempty(wslope)), wslope=500; end
if (nargin<3 || isempty(maxfe)),  maxfe=200; end

base = modpar26(1);
pv0  = getpar_l(base);
flst = logspace(log10(300), log10(6000), 160);   % moderate grid for slope (speed)
obj  = @(x) objf(x, pv0, fitidx, base, wslope, flst);

R0 = fdm26(struct('pa',base));
J0 = obj(pv0(fitidx));
fprintf('start: d=%.3f  maperr=%.1f  J=%.1f   (fitting %d params, wslope=%g, maxfe=%d)\n', ...
        R0.d, R0.maperr, J0, numel(fitidx), wslope, maxfe);

op = optimset('Display','off','MaxFunEvals',maxfe,'MaxIter',maxfe);
t0 = tic; [xopt,Jopt] = fminsearch(obj, pv0(fitidx), op); wall=toc(t0);

pv = pv0; pv(fitidx) = xopt; pa = setpar_l(base, pv);
Rf = fdm26(struct('pa',pa));                       % full-grid slope + maperr
evalc('S = tdm26(''coupeig'', struct(''pa'',pa));');
sub = S.maxRe_osc < 0;
fprintf('end:   d=%.3f  maperr=%.1f  J=%.1f   osc-maxRe=%+.1f (%s)   [%.0fs]\n', ...
        Rf.d, Rf.maperr, Jopt, S.maxRe_osc, tern(sub,'sub-critical','UNSTABLE'), wall);
fprintf('  d: %.3f -> %.3f (target 0.413) | maperr: %.1f -> %.1f (keep low) | stable: %d\n', ...
        R0.d, Rf.d, R0.maperr, Rf.maperr, sub);
% show which params moved
nm = parnames();
fprintf('  params moved: ');
for i = fitidx, fprintf('%s %.4g->%.4g  ', nm{i}, pv0(i), pv(i)); end
fprintf('\n');
R.pa=pa; R.pv=pv; R.pv0=pv0; R.R0=R0; R.Rf=Rf; R.S=S; R.fitidx=fitidx; R.wslope=wslope;
end

% ---- combined objective ----
function J = objf(x, pv0, fitidx, base, wslope, flst)
pv = pv0; pv(fitidx) = x; pa = setpar_l(base, pv);
try
    R = fdm26(struct('pa',pa,'flst',flst));
    J = R.maperr + wslope*abs(R.d - 0.413);
    if (~isfinite(J)), J = 1e9; end
catch
    J = 1e9;
end
end

% ---- 30-param vector <-> pa (getpar/setpar order, from fdm26/parfit24) ----
function nm = parnames()
nm = {'k1o','r1o','m1o','k2o','r2o','m2o','k3o','r3o','k4o','aco', ...
      'k1e','r1e','m1e','k2e','r2e','m2e','k3e','r3e','k4e','ace', ...
      'k1q','r1q','m1q','k2q','r2q','m2q','k3q','r3q','k4q','acq'};
end
function pv = getpar_l(pa),  nm=parnames(); pv=zeros(1,30); for i=1:30, pv(i)=pa.(nm{i}); end, end
function pa = setpar_l(pa,pv), nm=parnames(); for i=1:30, pa.(nm{i})=pv(i); end, end
function s = tern(c,a,b), if c, s=a; else, s=b; end, end
