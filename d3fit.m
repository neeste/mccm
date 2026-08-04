% D3FIT -- does the INTERNAL THIRD DOF pay for itself at nch=3?
%
% ===================== THE QUESTION, STATED PRECISELY ==================
%
% d3int=1 turns on the internal third DOF (no fluid compartment, no back-action;
% modpar26.m:14-27). From the fitted 2-DOF parameters of fit_nch3_surface it
% costs a lot of map accuracy, measured 2026-08-04 by d3pre:
%
%     d3int=0   maperr 144.02          <- the fitted 2-DOF frontier
%     d3int=1   maperr 948.76          <- the same parameters, third DOF live
%
% The third DOF is genuinely LIVE now -- k5o x2 moves maperr by +355 THROUGH THE
% PV VECTOR (not merely through pa), so the historical "0.000e+00 inertness" is
% gone and these six parameters are real coordinates. The question is whether a
% refit that is allowed to move them recovers the 144, and what it buys.
%
% ==================== WHY THERE ARE TWO ARMS, NOT ONE ==================
%
% This run adds r2e (15) and k3e (17) to the fit vector -- the two strongest
% latency levers, never in any default vector. So a one-armed run cannot
% attribute its result: if d3int=1 lands near 144, that could be the third DOF
% earning its keep, or it could be the two new levers, and the recorded history
% of this project is largely a history of exactly that confusion.
%
% So: TWO ARMS, SAME warm start, SAME objective, SAME vector -- differing ONLY
% in d3int and the six third-DOF coordinates that d3int makes meaningful.
%
%     ARM A   d3int=1   25 params  (18 + hbsc + k5o r5o m5o k5e r5e m5e)
%     ARM B   d3int=0   19 params  (18 + hbsc)                    <- CONTROL
%
% BUDGET-MATCHED IN EVALUATIONS, NOT HOURS. fminsearch makes progress per
% function evaluation; matching wall-clock would hand the cheaper arm 22% more
% search and make the comparison a statement about coupeig's cost. Measured by
% d3pre: 103 s/eval at d3int=1, 84 s at d3int=0 (coupeig 8.6 -> 23.9 s as the
% operator grows; abr_metric and score26 are near-flat).
%
%     ARM A   2 x 150 = 300 evals x 103 s = 8.6 h
%     ARM B   2 x 150 = 300 evals x  84 s = 7.0 h
%     total ~15.6 h + per-segment scoring, ~16.5 h inside a 20 h window
%
% The ~3.5 h margin is deliberate and not padding: two fits in this project have
% died silently to machine sleep, which is why every segment saves before the
% next begins. A lost segment costs 4.3 h, and the margin absorbs one.
%
% ARM A RUNS FIRST. It is the primary result, and a partial A with the control
% missing is worth more than a complete control with no A.
%
% ======================= WHAT COUNTS AS AN ANSWER ======================
%
% Read maperr and the shoulder together, in that order:
%
%   A reaches ~144 or better, B does not improve on it   -> the third DOF pays
%     for itself: it recovers what it cost AND buys the RL coordinate.
%   A and B both reach ~144                              -> r2e/k3e did the work.
%     The third DOF is affordable but not yet earning; the result is the two
%     new levers, and it should be reported as such.
%   A stalls well above B                                -> the third DOF costs
%     map accuracy it cannot recover at this budget. A real result, and the one
%     that says d3 needs its compartment (m=4) rather than the internal form.
%
% In every branch the SHOULDER COLUMN is the admissibility gate, not a nicety.
% If it rose across segments the run bought latency by growing the WNR's second
% peak and the latency numbers are not evidence -- see the header of s3fit.m for
% the measurement that established this.
% =======================================================================
RUNDIR = '/Users/neely/mccm_runs';
if (~exist(RUNDIR,'dir')), mkdir(RUNDIR); end
NSEG = 2; FE = 150;
WARM = 'fit_nch3_surface.mat';

% DIARY. matlab -batch BUFFERS stdout when redirected to a file and flushes only
% at exit, so a 16 h run is otherwise completely invisible until it is over --
% measured 2026-08-04, not assumed. diary writes incrementally, which is the
% difference between watching this run and hoping about it.
diary(fullfile(RUNDIR,'d3fit.log')); diary on;
cln = onCleanup(@() diary('off'));

HBMODE    = '';     % SHEAR drive, v1-v2 (SN's call, 2026-08-03). '' means "leave
                    % the model's own default", which is the shear branch --
                    % parfit26 only assigns pa.hbmode when hbm is non-empty.
                    % Asserted below rather than assumed, because a stale hbmode
                    % on the warm start would silently override it.
WSHOULDER = 2.0;    % as s3fit, derived there from the measured worth of the
                    % double-peak exploit (~0.09 of slope ~ 0.3-0.5 of Jsurf), so
                    % w>1.0 makes it unprofitable; 2.0 is a 2x margin. This run
                    % starts at shoulder 0.4103, INSIDE the region the guard has
                    % to pull out of, though less deeply than s3fit's 0.6195.
WSURF     = 1.0;
FITIDX    = [1 2 3 4 5 7 8 9 10 11 13 20 21 31 32 33 15 17];
                    % the s3fit vector PLUS r2e (15) and k3e (17). Both arms get
                    % them -- that is what makes the control a control.

L=load(WARM); pa0=L.R.pa;
if (isfield(pa0,'hbmode')), pa0=rmfield(pa0,'hbmode'); end

% ---- ASSERTIONS. This project's dominant failure mode is a run that reports
% ---- one configuration and fits another; assert instead of assuming.
if (isempty(HBMODE) && isfield(pa0,'hbmode') && ~isempty(pa0.hbmode))
    error('d3fit ABORTED: HBMODE is shear but the warm start carries hbmode="%s".', pa0.hbmode);
end
if (isfield(pa0,'latsoft') && isinf(pa0.latsoft))
    error('d3fit ABORTED: warm start sets latsoft=Inf (legacy gameable detector).');
end
if (isfield(pa0,'ohctau') && ~isempty(pa0.ohctau) && pa0.ohctau > 0)
    error(['d3fit ABORTED: warm start has the OHC RC pole ON (ohctau=%g). Measured ' ...
           '2026-08-04 to ABSORB energy above its corner (ohcBM sign-flips); it ' ...
           'must not be fitted in this form.'], pa0.ohctau);
end

fprintf('\n===== D3FIT: 2 arms x %d segments x %d evals =====\n', NSEG, FE);
fprintf('  A  d3int=1  25 params  ~%.1f h   (103 s/eval, measured)\n', NSEG*FE*103/3600);
fprintf('  B  d3int=0  19 params  ~%.1f h   ( 84 s/eval, measured)   CONTROL\n', NSEG*FE*84/3600);
fprintf('  objective: surface (16-cell band, free Delta) + shoulder guard + gain guard\n');
fprintf('  warm %s | wsurf %.2f | wshoulder %.2f | maptol 150 (nch=3 default)\n', ...
        WARM, WSURF, WSHOULDER);
fprintf('  output: %s\n', RUNDIR);

ARM = struct('tag',{'A','B'}, 'd3int',{1,0}, 'fitd3',{true,false});
res = struct('arm',{},'seg',{},'J',{},'slope',{},'lvlc',{},'sho',{}, ...
             'maperr',{},'amp',{},'osc',{},'hrs',{});
T0=tic;
for a=1:numel(ARM)
    PIN = struct('d3int', ARM(a).d3int);
    pa  = pa0; pa.d3int = ARM(a).d3int;

    % PRE-FLIGHT PER ARM. c3fit lost 7.3 h to a warm start parfit26 reported but
    % did not fit. The mapping differs between the arms -- setpar_l writes the six
    % third-DOF entries in both, but only arm A's model reads them -- so the check
    % is run per arm rather than once.
    H=parfit26('handles'); b0=modpar26(3); b0.d3int=ARM(a).d3int;
    pe = H.setpar(b0, H.getpar(pa), []);
    ew = fdm26(struct('pa',pa));  ew = ew.maperr;
    ee = fdm26(struct('pa',pe));  ee = ee.maperr;
    fprintf('\n===== ARM %s (d3int=%d) =====\n', ARM(a).tag, ARM(a).d3int);
    fprintf('  PRE-FLIGHT  warm %.4f | as jointobj builds it %.4f | delta %.3e\n', ew, ee, abs(ew-ee));
    if (abs(ew-ee) > 1)
        fprintf('  ARM %s SKIPPED: warm start does not survive setpar_l.\n', ARM(a).tag);
        continue;
    end

    for s=1:NSEG
        t0=tic; out=fullfile(RUNDIR,sprintf('d3fit_%s_seg%d.mat',ARM(a).tag,s));
        fprintf('----- arm %s segment %d/%d -----\n', ARM(a).tag, s, NSEG);
        try
            R = parfit26(3, struct('maxfe',FE, 'surface',1, 'wsurf',WSURF, ...
                                   'wshoulder',WSHOULDER, 'hbmode',HBMODE, ...
                                   'wgain',0.01, 'fithbsc',true, 'fitd3',ARM(a).fitd3, ...
                                   'fitidx',FITIDX, 'pin',PIN, 'warm',pa, ...
                                   'verbterm',true, 'out',out));
        catch e
            fprintf('  ARM %s SEGMENT %d FAILED: %s\n', ARM(a).tag, s, e.message(1:min(150,end)));
            break;
        end
        pa = R.pa;
        S = score26(pa,'fast',false);
        pm = pa; if (~isempty(HBMODE)), pm.hbmode = HBMODE; end
        m = abr_metric(pm,false);
        res(end+1)=struct('arm',ARM(a).tag,'seg',s,'J',R.J,'slope',m.slope, ...
            'lvlc',m.level_c,'sho',m.shoulder,'maperr',S.maperr,'amp',S.amp_gain, ...
            'osc',S.maxRe_osc,'hrs',toc(t0)/3600); %#ok<SAGROW>
        fprintf('\n  ARM %s SEG %d | J %.4f | slope %.4f | lvl_c %.3f | SHOULDER %.4f\n', ...
                ARM(a).tag, s, R.J, m.slope, m.level_c, m.shoulder);
        fprintf('              | maperr %.2f | amp %+.2f | osc %.1f | %.1f h\n', ...
                S.maperr, S.amp_gain, S.maxRe_osc, res(end).hrs);
        save(fullfile(RUNDIR,'d3fit_hist.mat'),'res');   % before the next segment:
    end                                                  % sleep costs one segment
end

fprintf('\n===== DONE, %.1f h total =====\n', toc(T0)/3600);
fprintf('  arm seg |      J | slope  | lvl_c | shldr  | maperr |   amp  |    osc | hrs\n');
for i=1:numel(res)
    h=res(i);
    fprintf('   %s   %d | %6.4f | %6.4f | %5.3f | %6.4f | %6.2f | %+6.2f | %6.1f | %4.1f\n', ...
            h.arm, h.seg, h.J, h.slope, h.lvlc, h.sho, h.maperr, h.amp, h.osc, h.hrs);
end
fprintf('\n  reference: d3int=0 unfitted 144.02 | d3int=1 unfitted 948.76\n');
fprintf('  READ THE SHOULDER COLUMN FIRST. If it rose, the latency numbers are not\n');
fprintf('  evidence. Then compare A vs B maperr: if BOTH reached ~144, the result is\n');
fprintf('  r2e/k3e, NOT the third DOF.\n');
fprintf('  saved: %s\n', fullfile(RUNDIR,'d3fit_hist.mat'));
disp('D3FIT_DONE');
