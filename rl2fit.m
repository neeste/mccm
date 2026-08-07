% RLFIT -- the LAST unexplored structural candidate: rlsplit, the RL carrying
% the fluid force instead of the TM.
%
% ============================ WHAT rlsplit IS ============================
%
% It does NOT add a fluid force to the reticular lamina. It MOVES one. Measured
% Df rows (rlprobe.m, 2026-08-05):
%
%   m=3 d3int    d1 [ 1  0 -1]   d2 [ 0 -1  1]   d3 ZERO   <- TM fluid-coupled
%   m=3 rlsplit  d1 [ 1 -1  0]   d2 ZERO         d3 [0 1 -1] <- RL fluid-coupled
%
% Both have exactly TWO fluid-coupled DOFs, because 3 chambers supply only 2
% independent pressure differences. The chambers become [ST, CL, SV], the BM
% spans ST-CL and the RL spans CL-SV: the two-membrane arrangement.
%
% ==================== THE RECORDED 517.10 IS FROM A BROKEN MODEL ==========
%
% "rlsplit: maperr 517.10 unfitted" has been carried as the baseline. It was
% measured with rlsplit but WITHOUT d3int, which is not a model: has3 is false,
% so ndof=2 while the rlsplit branch writes a third Df row. MATLAB grows Df to
% 3 rows, the solver's nd=min(...,ndof) ignores it, and the TM ends up clamped
% with nothing replacing it. fdm26 returns 517.10 anyway; tdm26 fails outright
% with a dimension error, so that configuration cannot time-march at all.
%
% THE CORRECT CONFIGURATION IS rlsplit + d3int, and it is well-posed:
%
%   maperr 338.00 | tdm wnr1 finite | maxRe +1.25 | osc -1080.9
%   slope 0.5945  | lvl_c 1.021 | shoulder 0.0626 (1/16 cells > 0.5)
%   amp   +8.21 dB
%
% Two things make this a better prospect than the internal third DOF was. It
% starts 194 points from the 2-DOF frontier of 144.02, not 805. And it starts
% nearly single-peaked, which cost d3int 150 evaluations to reach.
%
% ======================== WHY wgain IS 0.05, NOT 0.01 ====================
%
% The amplifier is the binding defect here, not the map. +8.21 dB against the
% project's stated bar of 40-60 dB is a model that has nearly stopped being a
% cochlea, and at the default weight the shortfall contributes 0.01*31.8 = 0.32
% to J, about a fifth of it, which is not enough to make it the priority. At
% 0.05 it contributes 1.59 and leads. DELIBERATE, and SN's call (2026-08-05):
% a rlsplit fit that lands a good map with a dead amplifier answers nothing.
% Every other weight is unchanged from d3fit so the runs stay comparable.
% =======================================================================
RUNDIR = '/Users/neely/mccm_runs';
if (~exist(RUNDIR,'dir')), mkdir(RUNDIR); end
diary(fullfile(RUNDIR,'rl2fit.log')); diary on;
cln = onCleanup(@() diary('off'));

NSEG = 2; FE = 150;
WSHOULDER = 2.0; WSURF = 1.0; HBMODE = ''; WGAIN = 0.05;
FITIDX = [1 2 3 4 5 7 8 9 10 11 13 20 21 31 32 33 15 17];   % as d3fit
PIN = struct('d3int',1, 'rlsplit',1, 'hbrl',-1);

L=load(fullfile(RUNDIR,'rlfit_seg2.mat')); pa=L.R.pa;   % the fitted rlsplit state
if (isfield(pa,'hbmode')), pa=rmfield(pa,'hbmode'); end
pa.d3int = 1; pa.rlsplit = 1; pa.hbrl = -1;

% ---- ASSERTIONS ----
if (isempty(HBMODE) && isfield(pa,'hbmode') && ~isempty(pa.hbmode))
    error('rlfit ABORTED: HBMODE is shear but the warm start carries hbmode="%s".', pa.hbmode);
end
if (isfield(pa,'latsoft') && isinf(pa.latsoft))
    error('rlfit ABORTED: warm start sets latsoft=Inf (legacy gameable detector).');
end
if (isfield(pa,'ohctau') && ~isempty(pa.ohctau) && pa.ohctau > 0)
    error('rlfit ABORTED: OHC RC pole is ON; measured 2026-08-04 to absorb above its corner.');
end
% THE ONE THAT MATTERS HERE. rlsplit without d3int is not a model, and it fails
% in the worst possible way: fdm26 returns a plausible maperr while tdm26 cannot
% run, so a fit would score the frequency domain of a model whose time domain
% does not exist. Assert rather than trust the pin struct.
if (~PIN.d3int)
    error('rlfit ABORTED: rlsplit requires d3int=1 (ndof would be 2 with a 3-row Df).');
end

fprintf('\n===== RL2FIT: rlsplit + d3int + hbrl=-1 (RL-driven MET), nch=3 =====\n');
fprintf('  %d segments x %d evals, 25 params, ~%.1f h at the measured 135 s/eval\n', ...
        NSEG, FE, NSEG*FE*135/3600);
fprintf('  start: maperr 338.00 | amp +8.21 dB | shoulder 0.0626 | slope 0.5945\n');
fprintf('  wgain %.2f (RAISED from 0.01: the amplifier is the binding defect)\n', WGAIN);
fprintf('  wsurf %.2f | wshoulder %.2f | maptol 150 | gainmin 40 dB\n', WSURF, WSHOULDER);
fprintf('  reference: 2-DOF frontier maperr 144.02 | d3int alone 948.76\n');

% ---- PRE-FLIGHT ----
H=parfit26('handles'); b0=modpar26(3);
pf=fieldnames(PIN); for i=1:numel(pf), b0.(pf{i})=PIN.(pf{i}); end
pe = H.setpar(b0, H.getpar(pa), []);
ew = fdm26(struct('pa',pa)); ew=ew.maperr;
ee = fdm26(struct('pa',pe)); ee=ee.maperr;
fprintf('\n  PRE-FLIGHT  warm %.4f | as jointobj builds it %.4f | delta %.3e\n', ew, ee, abs(ew-ee));
if (abs(ew-ee) > 1)
    error('rlfit ABORTED: warm start does not survive setpar_l (%.2f vs %.2f).', ew, ee);
end


% AMP GUARD. score26 returned amp_gain = NaN under hbrl=+1. wgain is the LEADING
% term in this objective by design, so a NaN there would make the largest term
% meaningless and the run worthless. Assert before spending 11 h.
s0 = score26(pa,'fast',false);
if (~isfield(s0,'amp_gain') || isempty(s0.amp_gain) || ~isfinite(s0.amp_gain(1)))
    error('rl2fit ABORTED: amp_gain is not finite at the start (%g).', ...
          subsref_d_(s0,'amp_gain'));
end
fprintf('  START amp %+.2f dB | maperr %.2f (sign settled: ohcBM +1.662e-04, injects)\n', ...
        s0.amp_gain(1), getfield(fdm26(struct('pa',pa)),'maperr'));

hist=struct('seg',{},'J',{},'slope',{},'lvlc',{},'sho',{},'maperr',{},'amp',{},'osc',{},'hrs',{});
T0=tic;
for s=1:NSEG
    t0=tic; out=fullfile(RUNDIR,sprintf('rl2fit_seg%d.mat',s));
    fprintf('----- segment %d/%d -----\n', s, NSEG);
    try
        R = parfit26(3, struct('maxfe',FE, 'surface',1, 'wsurf',WSURF, ...
                               'wshoulder',WSHOULDER, 'hbmode',HBMODE, ...
                               'wgain',WGAIN, 'fithbsc',true, 'fitd3',true, ...
                               'fitidx',FITIDX, 'pin',PIN, 'warm',pa, ...
                               'verbterm',true, 'out',out));
    catch e
        fprintf('  SEGMENT %d FAILED: %s\n', s, e.message(1:min(200,end))); break;
    end
    pa = R.pa; S = score26(pa,'fast',false);
    pm = pa; if (~isempty(HBMODE)), pm.hbmode = HBMODE; end
    m = abr_metric(pm,false);
    hist(end+1)=struct('seg',s,'J',R.J,'slope',m.slope,'lvlc',m.level_c, ...
        'sho',m.shoulder,'maperr',S.maperr,'amp',S.amp_gain,'osc',S.maxRe_osc, ...
        'hrs',toc(t0)/3600); %#ok<SAGROW>
    fprintf('\n  RL2FIT SEG %d | J %.4f | slope %.4f | lvl_c %.3f | SHOULDER %.4f\n', ...
            s, R.J, m.slope, m.level_c, m.shoulder);
    fprintf('              | maperr %.2f | AMP %+.2f | osc %.1f | %.1f h\n', ...
            S.maperr, S.amp_gain, S.maxRe_osc, hist(end).hrs);
    save(fullfile(RUNDIR,'rl2fit_hist.mat'),'hist');
end

fprintf('\n===== RL2FIT DONE, %.1f h =====\n', toc(T0)/3600);
fprintf('  seg |      J | slope  | lvl_c | shldr  | maperr |   amp  |    osc | hrs\n');
for i=1:numel(hist)
    h=hist(i);
    fprintf('  %3d | %6.4f | %6.4f | %5.3f | %6.4f | %6.2f | %+6.2f | %6.1f | %4.1f\n', ...
            h.seg, h.J, h.slope, h.lvlc, h.sho, h.maperr, h.amp, h.osc, h.hrs);
end
fprintf('\n  READ THE AMP COLUMN FIRST this time, not the shoulder. The model starts\n');
fprintf('  at +8.21 dB against a 40-60 dB bar; a good map with a dead amplifier is\n');
fprintf('  not a result. Then shoulder, then maperr against 144.02.\n');
disp('RL2FIT_DONE');
function v=subsref_d_(s,f), if (isfield(s,f)&&~isempty(s.(f))), v=double(s.(f)(1)); else, v=NaN; end, end
