% LCPROBE -- the slope-vs-level_c tradeoff at nch=3, measured, no fit.
%
% QUESTION. c3fit reached slope 0.4151 but walked level_c 4.796 -> 6.469
% (1.580 -> 1.885 %/dB against a 1.62 target) because wlevel=0 leaves level_c
% free anywhere inside the [3.5 6.5] band. SN keeps the band. So the question is
% no longer how to score level_c -- it is whether nch=3 CAN hold the slope while
% level_c comes back toward 5, or whether the two trade one-for-one the way ace
% traded them at nch=1 (see abr-tuning-levers: "Jslp and Jlcb trade nearly
% one-for-one ... the basin is FLAT by construction").
%
% WHY THIS IS CHEAP. Neither hbsc nor hbmx is visible to fdm26 (maperr) or to
% coupeig (maxRe_osc, which pins hbnl=0 itself at tdm26.m:1615). So the sweep
% needs ONLY abr_metric, ~115 s/point at nch=3, not the 184.7 s full evaluation.
% That invariance is ASSERTED below rather than assumed -- it is the whole basis
% of the cost estimate, and this project's dominant failure mode is exactly a
% second copy of a fact that agrees only while something is inert.
%
% ANCHOR. Row 1 reproduces c3fit_seg2 untouched and is checked against the
% values c3fit logged (slope 0.4151, level_c 6.469). If the anchor misses, the
% probe is not trustworthy and nothing below it means anything. (capstone-note
% failure #5: build the anchor into every decomposition.)
%
% SAFE TO RUN ALONGSIDE c3fit: tbabr.txt is written only when figures are
% requested (tdm26.m:938) and abr_metric(pa,false) requests none.
PROJ = '/Users/neely/Library/CloudStorage/OneDrive-FatherFlanagan''sBoysHome/STN/mod/mccm';
addpath(PROJ);
OUT  = '/Users/neely/mccm_runs/lcprobe.mat';

S2 = load('/Users/neely/mccm_runs/c3fit_seg2.mat');
pa = S2.R.pa;
fprintf('\n=== LCPROBE nch=3, from c3fit_seg2 ===\n');
fprintf('  m %d  d3int %d  hbsc %.6f  hbmx %.3e  ace %+.6f\n', ...
        pa.m, pa.d3int, pa.hbsc, pa.hbmx, pa.ace);

pct = @(lc) 100*(lc.^(1/100) - 1);   % level_c -> %/dB, target 1.62

% ---- ANCHOR ----
t0 = tic;
m0 = abr_metric(pa, false);
fprintf('\n  ANCHOR  slope %.4f (log 0.4151)  level_c %.3f (log 6.469)  %.3f %%/dB  [%.0f s]\n', ...
        m0.slope, m0.level_c, pct(m0.level_c), toc(t0));
if (abs(m0.slope-0.4151) > 0.002 || abs(m0.level_c-6.469) > 0.02)
    error('LCPROBE ANCHOR FAILED: seg2 does not reproduce its own logged values.');
end
fprintf('  ANCHOR OK.\n');

% ---- INVARIANCE ASSERTION: hbsc/hbmx must not reach fdm26 or coupeig ----
R0 = fdm26(struct('pa',pa));
evalc('C0 = tdm26(''coupeig'',struct(''pa'',pa));');
pt = pa; pt.hbsc = 0.25; pt.hbmx = 3e-9;
Rt = fdm26(struct('pa',pt));
evalc('Ct = tdm26(''coupeig'',struct(''pa'',pt));');
dmap = abs(R0.maperr - Rt.maperr);  dosc = abs(C0.maxRe_osc - Ct.maxRe_osc);
fprintf('  INVARIANCE  maperr %.4f vs %.4f (d %.2e) | osc %.3f vs %.3f (d %.2e)\n', ...
        R0.maperr, Rt.maperr, dmap, C0.maxRe_osc, Ct.maxRe_osc, dosc);
if (dmap > 1e-9 || dosc > 1e-9)
    error(['LCPROBE: hbsc/hbmx DO reach fdm26 or coupeig (dmap %.2e, dosc %.2e). ' ...
           'Every point below needs a full evaluation, not abr_metric alone.'], dmap, dosc);
end
fprintf('  INVARIANCE OK -- maperr %.2f and osc %.1f are CONSTANT over this sweep.\n', ...
        R0.maperr, C0.maxRe_osc);

% ---- SWEEPS ----
hb_sc = [0.005 0.010 0.020 0.030 pa.hbsc 0.060 0.100 0.200];
hb_mx = [1e-9 2e-9 4e-9 pa.hbmx 1e-8 2e-8 6e-8];
res = struct('lever',{},'val',{},'slope',{},'lvlc',{},'pdb',{},'rmse',{},'ok',{});

for i = 1:numel(hb_sc)
    p = pa; p.hbsc = hb_sc(i); t = tic;
    m = abr_metric(p, false);
    res(end+1) = struct('lever','hbsc','val',hb_sc(i),'slope',m.slope, ...
        'lvlc',m.level_c,'pdb',pct(m.level_c),'rmse',m.rmse,'ok',m.ok); %#ok<SAGROW>
    fprintf('  hbsc %8.5f | slope %.4f | lvl_c %7.3f | %5.3f %%/dB | ok %d | %.0f s\n', ...
            hb_sc(i), m.slope, m.level_c, pct(m.level_c), m.ok, toc(t));
    save(OUT,'res','pa','R0','C0');
end
for i = 1:numel(hb_mx)
    p = pa; p.hbmx = hb_mx(i); t = tic;
    m = abr_metric(p, false);
    res(end+1) = struct('lever','hbmx','val',hb_mx(i),'slope',m.slope, ...
        'lvlc',m.level_c,'pdb',pct(m.level_c),'rmse',m.rmse,'ok',m.ok); %#ok<SAGROW>
    fprintf('  hbmx %8.2e | slope %.4f | lvl_c %7.3f | %5.3f %%/dB | ok %d | %.0f s\n', ...
            hb_mx(i), m.slope, m.level_c, pct(m.level_c), m.ok, toc(t));
    save(OUT,'res','pa','R0','C0');
end

% ---- TABLE. J is the REAL parfit26 objective, band form, SN's choice. ----
fprintf('\n=== slope-vs-level_c at nch=3 (maperr %.2f, osc %.1f, amp fixed) ===\n', ...
        R0.maperr, C0.maxRe_osc);
fprintf('  lever      value | slope  | lvl_c  | %%/dB  | Jslp   Jlcb   |   J\n');
for i = 1:numel(res)
    r = res(i);
    Js = abs(r.slope - 0.413);
    Jl = 0.10*(max(0, 3.5 - r.lvlc) + max(0, r.lvlc - 6.5));
    fprintf('  %-5s %10.4g | %6.4f | %6.3f | %5.3f | %6.4f %6.4f | %6.4f\n', ...
            r.lever, r.val, r.slope, r.lvlc, r.pdb, Js, Jl, Js+Jl);
end
fprintf('\n  saved %s\n', OUT);
disp('LCPROBE_DONE');
