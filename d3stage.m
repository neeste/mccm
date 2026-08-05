% D3STAGE -- is the third DOF's 6 percent recovery a PHYSICS limit or a BUDGET limit?
%
% d3fit answered its question: turning on the internal third DOF costs 805 points
% of maperr and a 300-evaluation, 25-parameter refit recovered 51 of them (6%),
% while the budget-matched control did not move at all (0.46% in 155 evals),
% which eliminates r2e/k3e as an explanation.
%
% What d3fit did NOT establish is WHY the recovery was small, and there are two
% very different explanations:
%
%   PHYSICS   the internal third DOF genuinely cannot sit at the 2-DOF frontier;
%             no budget recovers the cost.
%   BUDGET    fminsearch in 25 dimensions from a point 2 J-units downhill is
%             simply too weak. The arms were matched in EVALUATIONS but not in
%             DISTANCE TO TRAVEL: arm B started at its own optimum and had
%             nothing to do, arm A started in a hole that d3int=1 had just dug.
%
% PHASE A settles this the cheap way, and it is the reason this script is not
% just "run arm A longer". Fit ONLY the six third-DOF coordinates: k5o r5o m5o
% k5e r5e m5e. Nothing else moves. A 6-dimensional simplex is far better
% conditioned than a 25-dimensional one (7 points to build instead of 26, and
% ~150 evaluations is a real budget in 6-D rather than a thin one in 25-D), and
% these six are the parameters that DOMINATE the cost: k5o alone moves maperr by
% 355. If six well-conditioned parameters cannot recover the 805 points, the
% limit is not the budget. If they recover a large share, then d3fit was simply
% under-budgeted and phase B exploits that.
%
% PHASE B then runs the full 25-parameter vector from whichever start is
% actually better, chosen by J rather than assumed: phase A's result or arm A's
% endpoint. Staging a badly-conditioned search behind a well-conditioned one is
% the standard remedy for exactly this failure mode.
%
% PHASE C is a probe, NOT a fit: where does the PROMOTED form (m=4, d3 with its
% own fluid compartment) actually sit? The Month-3 write-up says explicitly that
% this result does not settle the m=4 question, and the m=4 warm start is only
% now trustworthy, since clvent was being silently reset until today's fix.
% Cheap, and it says whether m=4 is worth a fit at all.
%
% COST, from measured rates (104 s/eval at d3int=1, arm A observed):
%   phase A   1 x 150 evals   ~4.3 h
%   phase B   2 x 150 evals   ~8.6 h
%   phase C   probe only      ~0.5 h
%   total                    ~13.4 h, inside a 20 h window with ~4 h margin
RUNDIR = '/Users/neely/mccm_runs';
if (~exist(RUNDIR,'dir')), mkdir(RUNDIR); end
diary(fullfile(RUNDIR,'d3stage.log')); diary on;
cln = onCleanup(@() diary('off'));

FE = 150;
WSHOULDER = 2.0; WSURF = 1.0; HBMODE = '';
FITIDX = [1 2 3 4 5 7 8 9 10 11 13 20 21 31 32 33 15 17];   % as d3fit
PIN = struct('d3int',1);
T0 = tic;

L=load('fit_nch3_surface.mat'); pa0=L.R.pa;
if (isfield(pa0,'hbmode')), pa0=rmfield(pa0,'hbmode'); end
pa0.d3int = 1;

fprintf('\n===== D3STAGE =====\n');
fprintf('  phase A: SIX third-DOF params only (well-conditioned)   %d evals\n', FE);
fprintf('  phase B: full 25-param vector from the better start   2 x %d evals\n', FE);
fprintf('  phase C: m=4 promoted-form probe (no fit)\n');

% ================= PHASE A: the six third-DOF coordinates alone =================
fprintf('\n===== PHASE A: k5o r5o m5o k5e r5e m5e ONLY =====\n');
RA = [];
try
    RA = parfit26(3, struct('maxfe',FE, 'surface',1, 'wsurf',WSURF, ...
                            'wshoulder',WSHOULDER, 'hbmode',HBMODE, 'wgain',0.01, ...
                            'fitidx',[], 'fitd3',true, 'fithbsc',false, ...
                            'pin',PIN, 'warm',pa0, 'verbterm',true, ...
                            'out',fullfile(RUNDIR,'d3stage_A.mat')));
    SA=score26(RA.pa,'fast',false); mA=abr_metric(RA.pa,false);
    fprintf('\n  PHASE A | J %.4f | maperr %.2f | shoulder %.4f | slope %.4f | lvl_c %.3f\n', ...
            RA.J, SA.maperr, mA.shoulder, mA.slope, mA.level_c);
    fprintf('  recovery: 948.76 -> %.2f, that is %.1f%% of the 805-point cost\n', ...
            SA.maperr, 100*(948.76-SA.maperr)/(948.76-144.02));
catch e
    fprintf('  PHASE A FAILED: %s\n', e.message(1:min(200,end)));
end

% ================= PHASE B: full vector from the better start ==================
fprintf('\n===== PHASE B: full 25-param vector =====\n');
% CHOOSE BY J, DO NOT ASSUME. Phase A optimizes the same objective, so its J is
% directly comparable to arm A's. Picking the wrong start would waste 8.6 h
% re-deriving something already in hand.
startpa = []; startJ = Inf; startwho = '';
fA = fullfile(RUNDIR,'d3fit_A_seg2.mat');
if (exist(fA,'file')), LA=load(fA); startpa=LA.R.pa; startJ=LA.R.J; startwho='d3fit arm A seg2'; end
if (~isempty(RA) && RA.J < startJ), startpa=RA.pa; startJ=RA.J; startwho='phase A'; end
if (isempty(startpa)), startpa=pa0; startJ=NaN; startwho='unfitted d3int=1'; end
fprintf('  start: %s (J %.4f)\n', startwho, startJ);
pa = startpa; if (isfield(pa,'hbmode')), pa=rmfield(pa,'hbmode'); end
pa.d3int = 1;
for s=1:2
    t0=tic;
    fprintf('----- phase B segment %d/2 -----\n', s);
    try
        RB = parfit26(3, struct('maxfe',FE, 'surface',1, 'wsurf',WSURF, ...
                                'wshoulder',WSHOULDER, 'hbmode',HBMODE, 'wgain',0.01, ...
                                'fitidx',FITIDX, 'fitd3',true, 'fithbsc',true, ...
                                'pin',PIN, 'warm',pa, 'verbterm',true, ...
                                'out',fullfile(RUNDIR,sprintf('d3stage_B_seg%d.mat',s))));
    catch e
        fprintf('  PHASE B SEG %d FAILED: %s\n', s, e.message(1:min(200,end))); break;
    end
    pa = RB.pa; S=score26(pa,'fast',false); m=abr_metric(pa,false);
    fprintf('\n  PHASE B SEG %d | J %.4f | maperr %.2f | SHOULDER %.4f | slope %.4f | lvl_c %.3f | %.1f h\n', ...
            s, RB.J, S.maperr, m.shoulder, m.slope, m.level_c, toc(t0)/3600);
    fprintf('               | amp %+.2f | maxRe_osc %.1f\n', S.amp_gain, S.maxRe_osc);
end

% ================= PHASE C: where does the PROMOTED form sit? ==================
% Probe only. The m=4 warm start became trustworthy only today: clvent was
% silently reset to its default on every evaluation until the mapping fix, so
% any earlier m=4 number was measured on a different model than the one saved.
fprintf('\n===== PHASE C: m=4 promoted-form probe (no fit) =====\n');
try
    L4=load('refit_m4_full.mat'); p4=L4.R.pa;
    if (isfield(p4,'hbmode')), p4=rmfield(p4,'hbmode'); end
    S4=score26(p4,'fast',false); m4=abr_metric(p4,false);
    fprintf('  m=4 at its saved fit (clvent %.3g, now correctly carried):\n', p4.clvent);
    fprintf('    maperr %.2f | amp %+.2f | maxRe %+.2f | osc %.1f\n', ...
            S4.maperr, S4.amp_gain, S4.maxRe, S4.maxRe_osc);
    fprintf('    slope %.4f | lvl_c %.3f | shoulder %.4f | %d/16 cells>0.5\n', ...
            m4.slope, m4.level_c, m4.shoulder, sum(m4.sho(isfinite(m4.sho))>0.5));
    fprintf('  compare: nch=3 2-DOF maperr 144.02 slope 0.4927 lvl_c 5.602 shoulder 0.2603\n');
catch e
    fprintf('  PHASE C FAILED: %s\n', e.message(1:min(200,end)));
end

fprintf('\n===== D3STAGE DONE, %.1f h =====\n', toc(T0)/3600);
disp('D3STAGE_DONE');
