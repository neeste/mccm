% S3FIT -- nch=3 against the 16-CELL LATENCY SURFACE, with the peak-switch guards on.
%
% ===================== WHY THE SCALAR SLOPE IS DEAD =====================
%
% c3fit reached slope 0.4118 against a 0.413 target and J 0.0012. That number is
% an artifact of the WNR peak detector, measured 2026-08-03 (drvcheck.m):
%
%     drive       slope    lvl_c   shoulder
%     default    0.4118    6.399     0.6195     <- what c3fit reported
%     bm         0.5058    6.244     0.4036
%     dof2       0.4861    6.808     0.5407
%
% |slope(bm) - slope(default)| = 0.0940, which is 78x the 0.0012 gap to target.
% The result is a statement about the latency drive, not about the cochlea.
%
% CAUSE. wnr_latency takes the LARGEST local maximum. When the WNR carries two
% comparable peaks the detector flips between them, and the reported latency
% jumps discontinuously in the parameters. At c3fit's optimum the shoulder ratio
% (2nd peak / main) is 0.6195 on average with 10/16 cells above 0.5 and three
% cells at 0.999-1.000 -- two peaks of EQUAL height, i.e. a coin flip. Compare
% the 60 dB column, default vs bm (ms):
%     default   4.92  3.18  2.22  2.20
%     bm       10.04  6.80  4.56  2.52
% The default drive reports roughly HALF the latency at three of four
% frequencies. That is the switch, and it is almost certainly the recorded
% "kink at 60 dB" in the latency surface, which has been carried as a physical
% finding.
%
% THE FIT CAUSES IT. Shoulder rises monotonically with fitting effort:
%     stock nch=1 0.0000 -> nch=1 longfit 0.3579
%     stock nch=3 0.3868 -> warm 0.4915 -> c3fit seg1 0.6147 -> seg3 0.6195
% So the optimizer is GROWING the second peak until it overtakes the first,
% because that shortens the reported latency and bends the slope to target.
% wnr_latency's own comment records this exploit once already, under a
% first-local-max rule; taking the largest peak did not remove it, it inverted
% it. Every fitted slope in this project is contaminated, and the contamination
% grows with the budget spent.
%
% stock nch=1 has shoulder 0.0000 over all 16 cells, so a single-peaked WNR is
% reachable -- this is a defect the objective can be asked to avoid, not a fact
% of the model.
%
% ========================= WHAT THIS RUN DOES ==========================
%
% 1. SURFACE, not slope. opts.surface=1 scores the raw 16-cell latency matrix
%    against the band spanning the 1988 and 2013 power laws, with the excess
%    neural delay Delta free in [0 1.5] ms (parfit26's surf_term). Uses all 16
%    cells instead of compressing them to two summary numbers, and does not
%    assume the surface is power-law shaped -- which the recorded 13% shape
%    residual says it is not.
%
%    NOTE surface=1 DROPS the level_c band term: wlcb lives in the other branch
%    of jointobj. level_c is then constrained only implicitly, through the level
%    axis of the surface. It is reported every segment so the drift that c3fit
%    hid (4.796 -> 6.399) stays visible.
%
% 2. wshoulder > 0. Without it the surface is scored from the same switching
%    detector and the same exploit is available -- the surface reads those 16
%    latencies. This is the guard that makes the objective well-posed, and it is
%    why this run is not simply c3fit with a different ABR term.
%
% 3. hbmode PINNED. Left as a knob because it is a physics choice, not a fix:
%    the default v1-v2 is the shear that actually drives the hair bundle, while
%    'bm' (v1 alone) is labelled a DIAGNOSTIC in tdm26 and discards the
%    second-DOF contribution the multi-chamber model exists to represent. 'bm'
%    still shows shoulder 0.4036, so it reduces the ambiguity without removing
%    it. SET DELIBERATELY -- see the note where HBMODE is defined below.
%
% Warm start c3fit_seg3.mat, pinned d3int=0, and the same setpar_l pre-flight
% c3fit added after 7.3 h were lost to a warm start parfit26 reported but did
% not fit.
%
% COST. abr_metric is now ~6.7 s (tbabr parfor), so an evaluation is dominated
% by score26 at ~69 s: roughly 76 s, against 184.7 s before. 3 x 128 = 384
% evaluations is about 8.1 h rather than 19.7 h.
% =======================================================================
RUNDIR = '/Users/neely/mccm_runs';
if (~exist(RUNDIR,'dir')), mkdir(RUNDIR); end
NSEG = 3; FE = 128;
PIN  = struct('d3int',0);

% ---- THE THREE DELIBERATE CHOICES ----
HBMODE    = '';     % SHEAR drive, v1-v2 (SN's call, 2026-08-03). This is what
                    % actually deflects the hair bundle, and it is the drive the
                    % multi-chamber model exists to represent -- 'bm' (v1 alone)
                    % is labelled a DIAGNOSTIC in tdm26 and throws away the
                    % second-DOF contribution. parfit26.m:183 recommends 'bm',
                    % but that recommendation was a workaround for the detector
                    % defect fixed in b340548; with a continuous detector the
                    % reason for it is gone. COST: the warm start sits at
                    % shoulder 0.6195 under shear vs 0.4036 under bm, so this
                    % run STARTS inside the double-peaked region and WSHOULDER
                    % has to pull it out rather than merely keep it out.
                    %
                    % NOTE parfit26 only assigns pa.hbmode when hbm is non-empty,
                    % so '' means "leave the model's own default", which is the
                    % shear branch in hc_step. Asserted below rather than assumed.
WSHOULDER = 2.0;    % weight on the WNR shoulder ratio, DERIVED not guessed.
                    % The exploit was worth ~0.09 of slope, which is roughly
                    % 0.3-0.5 ms of band violation, i.e. 0.3-0.5 of Jsurf. For
                    % double-peaking to be unprofitable the shoulder term must
                    % cost more than that over the range it moves (~0.5), so
                    % w > 1.0. 2.0 gives a 2x margin over the measured value of
                    % the exploit. Start scale (calibration, shear drive):
                    % surf_rms 0.609 ms, shoulder 0.6195 -> Jsurf 0.61 vs
                    % Jsho 1.24, so the guard leads early and should recede as
                    % the shoulder comes down.
WSURF     = 1.0;

WARM = fullfile(RUNDIR,'c3fit_seg3.mat');   % outputs live in RUNDIR, not the repo
L=load(WARM); if (isfield(L,'R')), pa=L.R.pa; else, pa=L.pa; end
pf=fieldnames(PIN); for i=1:numel(pf), pa.(pf{i})=PIN.(pf{i}); end

% ---- PRE-FLIGHT: does the warm start survive the parameter mapping? ----
H=parfit26('handles'); b0=modpar26(3);
for i=1:numel(pf), b0.(pf{i})=PIN.(pf{i}); end
pa_eval = H.setpar(b0, H.getpar(pa), []);
e_warm = fdm26(struct('pa',pa));      e_warm = e_warm.maperr;
e_eval = fdm26(struct('pa',pa_eval)); e_eval = e_eval.maperr;
fprintf('\n  PRE-FLIGHT  warm-start maperr %.2f | as jointobj builds it %.2f | delta %.2f\n', ...
        e_warm, e_eval, abs(e_warm-e_eval));
if (abs(e_warm - e_eval) > 1)
    error(['s3fit ABORTED: the warm start does not survive setpar_l ' ...
           '(maperr %.2f vs %.2f).'], e_warm, e_eval);
end
fprintf('  PRE-FLIGHT PASSED.\n');

% ---- DRIVE / DETECTOR ASSERTIONS ----
% parfit26 leaves pa.hbmode alone when hbm is '', so a STALE hbmode carried on
% the warm start would silently override the choice made above and the run would
% report one drive while fitting another. That is this project's dominant
% failure mode; assert instead of assuming.
if (isempty(HBMODE) && isfield(pa,'hbmode') && ~isempty(pa.hbmode))
    error('s3fit ABORTED: HBMODE is shear but the warm start carries hbmode="%s".', pa.hbmode);
end
% latsoft Inf would restore the gameable largest-peak detector, which is the
% whole reason this run exists.
if (isfield(pa,'latsoft') && isinf(pa.latsoft))
    error('s3fit ABORTED: warm start sets latsoft=Inf (legacy gameable detector).');
end
fprintf('  DRIVE: %s | detector: soft-argmax softp=%g\n', ...
        ternstr_l(isempty(HBMODE),'shear (v1-v2)',HBMODE), ...
        subsref_default_l(pa,'latsoft',8));

% ---- START BREAKDOWN, both drives, so the weights are set from measurement ----
pd = pa; if (isfield(pd,'hbmode')), pd=rmfield(pd,'hbmode'); end
pb = pa; pb.hbmode = 'bm';
m0 = abr_metric(pd,false); m1 = abr_metric(pb,false);
fprintf('\n  START  drive      slope   lvl_c   shoulder   cells>0.5\n');
fprintf('         default   %7.4f %7.3f   %8.4f   %2d/16\n', ...
        m0.slope, m0.level_c, m0.shoulder, sum(m0.sho(isfinite(m0.sho))>0.5));
fprintf('         bm        %7.4f %7.3f   %8.4f   %2d/16\n', ...
        m1.slope, m1.level_c, m1.shoulder, sum(m1.sho(isfinite(m1.sho))>0.5));

fprintf('\n===== S3FIT nch=3: %d segments x %d evals (~%.1f h at 76 s/eval) =====\n', ...
        NSEG, FE, NSEG*FE*76/3600);
fprintf('  objective: surface (16-cell band, free Delta) + shoulder guard + gain guard\n');
fprintf('  hbmode "%s" | wsurf %.2f | wshoulder %.2f | maptol 150 | wgain 0.01\n', ...
        HBMODE, WSURF, WSHOULDER);
fprintf('  warm: c3fit_seg3.mat  (shear drive, FIXED detector: slope 0.4409,\n');
fprintf('        level_c 7.056 -- OUT of the old [3.5 6.5] band; it read 0.4118\n');
fprintf('        and 6.399 through the gameable detector. shoulder 0.6195)\n');
fprintf('  output: %s\n\n', RUNDIR);

hist=struct('seg',{},'J',{},'slope',{},'lvlc',{},'sho',{},'maperr',{},'amp',{},'osc',{},'hrs',{});
bestJ=Inf; bestpa=pa; T0=tic;
for s=1:NSEG
    t0=tic; out=fullfile(RUNDIR,sprintf('s3fit_seg%d.mat',s));
    fprintf('----- segment %d/%d -----\n', s, NSEG);
    try
        R = parfit26(3, struct('maxfe',FE, 'surface',1, 'wsurf',WSURF, ...
                               'wshoulder',WSHOULDER, 'hbmode',HBMODE, ...
                               'wgain',0.01, 'fithbsc',true, 'pin',PIN, ...
                               'warm',pa, 'verbterm',true, 'out',out));
    catch e
        fprintf('  SEGMENT %d FAILED: %s\n', s, e.message(1:min(150,end)));
        break;
    end
    pa = R.pa;
    S = score26(pa,'fast',false);
    pm = pa; if (~isempty(HBMODE)), pm.hbmode = HBMODE; end
    m = abr_metric(pm,false);
    hrs=toc(t0)/3600;
    J = R.J;
    if (J<bestJ), bestJ=J; bestpa=pa; end
    hist(end+1)=struct('seg',s,'J',J,'slope',m.slope,'lvlc',m.level_c, ...
        'sho',m.shoulder,'maperr',S.maperr,'amp',S.amp_gain,'osc',S.maxRe_osc, ...
        'hrs',hrs); %#ok<SAGROW>
    fprintf('\n  SEG %d | J %.4f | slope %.4f | lvl_c %.3f | SHOULDER %.4f\n', ...
            s, J, m.slope, m.level_c, m.shoulder);
    fprintf('        | maperr %.2f | amp %+.2f | osc %.1f | %.1f h\n', ...
            S.maperr, S.amp_gain, S.maxRe_osc, hrs);
    save(fullfile(RUNDIR,'s3fit_hist.mat'),'hist','bestJ','bestpa');
end

fprintf('\n===== DONE, %.1f h total =====\n', toc(T0)/3600);
fprintf('  seg |      J | slope  | lvl_c | shldr  | maperr |   amp  |    osc | hrs\n');
for i=1:numel(hist)
    h=hist(i);
    fprintf('  %3d | %6.4f | %6.4f | %5.3f | %6.4f | %6.2f | %+6.2f | %6.1f | %4.1f\n', ...
            h.seg, h.J, h.slope, h.lvlc, h.sho, h.maperr, h.amp, h.osc, h.hrs);
end
fprintf('\n  READ THE SHOULDER COLUMN FIRST. If it rose across segments, the fit is\n');
fprintf('  still buying latency with the second peak and the slope is not evidence.\n');
fprintf('  saved: %s\n', fullfile(RUNDIR,'s3fit_hist.mat'));
disp('S3FIT_DONE');

% local helpers (script-local; MATLAB requires them at the end of the file)
function s=ternstr_l(c,a,b), if c, s=a; else, s=b; end, end
function v=subsref_default_l(s,f,d), if (isfield(s,f)), v=s.(f); else, v=d; end, end
