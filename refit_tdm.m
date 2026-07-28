function R = refit_tdm(nch, fitidx, maxfe, outfile, warmfile, wmap, lclo, lchi)
% REFIT_TDM  tdm-ANCHORED joint re-fit (1- or 3-chamber). Objective = true
%   time-domain ABR slope, with tuning, sub-criticality, and level-dependence
%   enforced in-loop:
%     J = |tdm_slope-0.413| + 0.005*[osc_maxRe+40]^+ + 0.001*[maperr-tol]^+ + 0.20*[2-level_c]^+
%   (tol = good-tuning target ~1.6x the 1-chamber baseline maperr, so a 3-chamber
%   fit is rewarded for reaching MAP AND forward latency together.)
%
%   nch     : 1 or 3 (selects par_CEL16 / modpar26c3).
%   fitidx  : indices into the 30-param impedance vector to fit.
%   maxfe   : fminsearch MaxFunEvals.
%   outfile : .mat to save R into.
%   warmfile: optional .mat whose R.pa warm-starts the fit (else baseline).

if (nargin<1 || isempty(nch)),     nch=1; end
if (nargin<2 || isempty(fitidx)),  fitidx=[1 2 3 11 12 13]; end
if (nargin<3 || isempty(maxfe)),   maxfe=60; end
if (nargin<4 || isempty(outfile)), outfile=sprintf('refit_tdm_nch%d.mat',nch); end
if (nargin<5), warmfile=''; end
if (nargin<6 || isempty(wmap)), wmap=0.0004; end
if (nargin<7 || isempty(lclo)), lclo=2.0; end
if (nargin<8 || isempty(lchi)), lchi=8.0; end

base = modpar26(nch);
pa0  = base;
if (~isempty(warmfile) && exist(warmfile,'file')), L=load(warmfile); pa0=L.R.pa; end
pv0  = getpar_l(pa0);

tol.map  = 1.6 * fdm26(struct('pa',modpar26(1))).maperr;   % good-tuning target (~167)
tol.crit = -40; tol.lclo = lclo; tol.lchi = lchi;          % level_c band (tighten to center level dependence)
tol.wmap = wmap;                                            % tuning-penalty weight (raise for a tuning-focused pass)
obj = @(x) objf_tdm(x, pv0, fitidx, base, tol);

m0 = abr_metric(base, false);
fprintf('nch=%d baseline: tdm_d=%.3f  level_c=%.2f  maperr=%.1f (tol %.0f)\n', nch, m0.slope, m0.level_c, fdm26(struct('pa',base)).maperr, tol.map);
J0 = obj(pv0(fitidx));
fprintf('start: J=%.4f  (fitting %d params, maxfe=%d, out=%s)\n', J0, numel(fitidx), maxfe, outfile);

op = optimset('Display','off','MaxFunEvals',maxfe,'MaxIter',maxfe,'OutputFcn',@outfun);
t0 = tic; [xopt,Jopt] = fminsearch(obj, pv0(fitidx), op); wall = toc(t0);

pv = pv0; pv(fitidx) = xopt; pa = setpar_l(base, pv);
mf = abr_metric(pa, false); evalc('S = tdm26(''coupeig'', struct(''pa'',pa));'); Rf = fdm26(struct('pa',pa));
fprintf('\n=== nch=%d tdm-anchored re-fit (%.0f s, J %.4f->%.4f) ===\n', nch, wall, J0, Jopt);
fprintf('  ABR slope d : %.3f -> %.3f   (target 0.413)\n', m0.slope, mf.slope);
fprintf('  level_c     : %.2f -> %.2f   (target 5)\n', m0.level_c, mf.level_c);
fprintf('  maperr      : %.1f -> %.1f   (good-tuning tol %.0f)\n', fdm26(struct('pa',base)).maperr, Rf.maperr, tol.map);
fprintf('  osc-maxRe   : %+.1f (%s), n_sub=%d\n', S.maxRe_osc, tern(S.maxRe_osc<0,'sub-critical','UNSTABLE'), mf.n_sub);
R.pa=pa; R.pv=pv; R.m0=m0; R.mf=mf; R.S=S; R.Rf=Rf; R.fitidx=fitidx; R.nch=nch;
save(outfile,'R');
end

function J = objf_tdm(x, pv0, fitidx, base, tol)
pv = pv0; pv(fitidx) = x; pa = setpar_l(base, pv);
try
    m = abr_metric(pa, false);
    if (~m.ok), J = 1e6; return; end
    evalc('S = tdm26(''coupeig'', struct(''pa'',pa));');
    Rf = fdm26(struct('pa',pa));
    J = abs(m.slope - 0.413) ...
      + 0.005  * max(0, S.maxRe_osc - tol.crit) ...
      + tol.wmap * max(0, Rf.maperr  - tol.map) ...
      + 0.10   * (max(0, tol.lclo - m.level_c) + max(0, m.level_c - tol.lchi));
    if (~isfinite(J)), J = 1e6; end
catch
    J = 1e6;
end
end

function stop = outfun(~, ov, state)
stop = false;
if (strcmp(state,'iter')), fprintf('  iter %2d: J=%.4f  (fevals %d)\n', ov.iteration, ov.fval, ov.funccount); end
end

function nm = parnames()
nm = {'k1o','r1o','m1o','k2o','r2o','m2o','k3o','r3o','k4o','aco', ...
      'k1e','r1e','m1e','k2e','r2e','m2e','k3e','r3e','k4e','ace', ...
      'k1q','r1q','m1q','k2q','r2q','m2q','k3q','r3q','k4q','acq'};
end
% params 1..30 = impedance vector; 31..(30+numel(chsz)) = chamber sizes (the 3-chamber structural lever)
function pv = getpar_l(pa)
nm=parnames(); nc=numel(pa.chsz); pv=zeros(1,30+nc);
for i=1:30, pv(i)=pa.(nm{i}); end
pv(31:30+nc)=pa.chsz(:)';
end
function pa = setpar_l(pa,pv)
nm=parnames(); nc=numel(pa.chsz);
for i=1:30, pa.(nm{i})=pv(i); end
pa.chsz=pv(31:30+nc);
end
function s = tern(c,a,b), if c, s=a; else, s=b; end, end
