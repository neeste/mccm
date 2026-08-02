% LONGFIT -- the Aug 3-9 run. nch=1, 16 params, ABR solved, CF MAP made active.
%
% ============================ WHY THIS SHAPE ============================
%
% THE ABR OBJECTIVE IS ALREADY SOLVED, so this run is NOT chasing it. At
% ace=+0.40 with hbsc=0.64 (both reachable from sweep_nch1b):
%
%     slope   0.4160   target 0.413        Jslp  0.0030
%     level_c 3.732    band [3.5 6.5]      Jlcb  0.0000   <- in band, a project first
%     maperr  90.52    maptol 105          Jmap  0.0000
%     amp     +44.62   gainmin 40          Jgain 0.0000
%     osc     -83.5    margin -40          Josc  0.0000
%                                          J     0.0030   (fitted optimum: 0.1928)
%
% Pointing 900 evaluations at a J of 0.0030 would mostly wander. SN's call
% (2026-08-01) was to make the CF MAP the active constraint instead: hold the ABR
% terms at zero and let the fit improve the tonotopic map underneath them.
%
% maptol 80, NOT the nch<=2 default of 105. At maperr 90.52 that puts
% Jmap = 0.001*10.52 = 0.0105, three times the slope term, so the map actually
% pulls -- and because the term is a HINGE, all pressure returns to the ABR the
% moment maperr reaches 80. That makes the goal bounded and well posed: get the
% map under 80 WITHOUT breaking slope or level_c. For scale, the best maperr ever
% measured in this project is 79.17 (ace=-0.10, acecheck) and the standing record
% before today was 90.97.
%
% The trade this permits is deliberate and roughly even: 0.01 of slope costs the
% same J as 10 of maperr. If the fit spends slope 0.416 -> 0.426 to buy maperr
% 90 -> 80, that is the intended behaviour, not a failure.
%
% 16 PARAMETERS: the 15 of the default nch=1 fitidx (which ALREADY contains ace
% at index 20 -- it was there all along and the 150-eval run moved it 0.20%),
% plus hbsc via opts.fithbsc. gampro is deliberately NOT fitted: measured
% CF-map-free but stability-bankrupt, ~2600 dB of oscillatory margin per unit of
% slope, so it costs a dimension and buys nothing.
%
% RESTARTS, NOT ONE LONG CALL. The 150-eval run flatlined at J=0.1928 from
% iteration 74 to 84 with ace essentially untouched: Nelder-Mead simplex collapse
% in 15 dimensions. Pouring 900 evaluations into one collapsed simplex reproduces
% that at four times the cost. Each segment re-inflates around the previous best.
%
% OUTPUT IS LOCAL, NOT ONEDRIVE (SN). parfit26 checkpoints EVERY iteration to
% [out '.ckpt.mat'], ~500 writes here; under OneDrive that is 500 upload cycles
% from a conference connection, with a sync conflict able to land on the file the
% run is actively writing.
%
% SLEEP. The laptop sleeps in transit. macOS suspends and resumes MATLAB cleanly,
% so the run survives, but the hours printed below are WALL clock and will read
% long by however long the lid was shut. Judge progress by evaluations, not hours.
%
% =======================================================================
RUNDIR = '/Users/neely/mccm_runs';
if (~exist(RUNDIR,'dir')), mkdir(RUNDIR); end
NSEG = 4; FE = 225;          % ~900 evals at 71.4 s/eval ~ 17.9 h
ACE0 = 0.40; HBSC0 = 0.64; MAPTOL = 80;

L=load('sweep_nch1b.mat'); pa=L.R.pa; pa.ace=ACE0; pa.hbsc=HBSC0;
fprintf('\n===== LONGFIT nch=1: %d segments x %d evals (~%.1f h) =====\n', ...
        NSEG, FE, NSEG*FE*71.4/3600);
fprintf('  start: ace %+.3f  hbsc %.3f   (J 0.0030; ABR solved)\n', ACE0, HBSC0);
fprintf('  maptol %d (not the default 105) -- the CF map is the active constraint\n', MAPTOL);
fprintf('  16 params: default nch=1 fitidx (ace is index 20) + hbsc\n');
fprintf('  output: %s\n\n', RUNDIR);

hist=struct('seg',{},'J',{},'slope',{},'lvlc',{},'maperr',{},'amp',{},'osc',{},'ace',{},'hbsc',{},'hrs',{});
bestJ=Inf; bestpa=pa; T0=tic;
for s=1:NSEG
    t0=tic; out=fullfile(RUNDIR,sprintf('longfit_seg%d.mat',s));
    fprintf('----- segment %d/%d -----\n', s, NSEG);
    try
        R = parfit26(1, struct('maxfe',FE, 'wgain',0.01, 'fithbsc',true, ...
                               'maptol',MAPTOL, 'warm',pa, 'verbterm',true, ...
                               'out',out));
    catch e
        fprintf('  SEGMENT %d FAILED: %s\n', s, e.message(1:min(150,end)));
        break;                       % every completed segment stays on disk
    end
    pa = R.pa;                       % chain: next simplex is built around this
    S = score26(pa,'fast',false); m = abr_metric(pa,false);
    Js=abs(m.slope-0.413);
    Jl=0.10*(max(0,3.5-m.level_c)+max(0,m.level_c-6.5));
    Jp=0.001*max(0,S.maperr-MAPTOL); Jg=0.01*max(0,40-S.amp_gain);
    Jo=0.005*max(0,S.maxRe_osc+40);
    J=Js+Jl+Jp+Jg+Jo; hrs=toc(t0)/3600;
    if (J<bestJ), bestJ=J; bestpa=pa; end
    hist(end+1)=struct('seg',s,'J',J,'slope',m.slope,'lvlc',m.level_c, ...
        'maperr',S.maperr,'amp',S.amp_gain,'osc',S.maxRe_osc, ...
        'ace',pa.ace,'hbsc',pa.hbsc,'hrs',hrs); %#ok<SAGROW>
    fprintf('\n  SEG %d | J %.4f | slope %.4f | lvl_c %.3f | maperr %.2f | amp %+.2f | osc %.1f\n', ...
            s, J, m.slope, m.level_c, S.maperr, S.amp_gain, S.maxRe_osc);
    fprintf('        | ace %+.4f | hbsc %.4f | %.1f h\n', pa.ace, pa.hbsc, hrs);
    fprintf('    terms: slope %.4f  lcb %.4f  map %.4f  gain %.4f  osc %.4f\n\n', Js, Jl, Jp, Jg, Jo);
    save(fullfile(RUNDIR,'longfit_hist.mat'),'hist','bestJ','bestpa');
end

fprintf('\n===== DONE, %.1f h total =====\n', toc(T0)/3600);
fprintf('  seg |      J | slope  | lvl_c | maperr |   amp  |    osc |    ace |  hbsc  | hrs\n');
for i=1:numel(hist)
    h=hist(i);
    fprintf('  %3d | %6.4f | %6.4f | %5.3f | %6.2f | %+6.2f | %6.1f | %+6.4f | %6.4f | %4.1f\n', ...
            h.seg, h.J, h.slope, h.lvlc, h.maperr, h.amp, h.osc, h.ace, h.hbsc, h.hrs);
end
fprintf('\n  best J %.4f    start 0.0030 (at maptol 105); 150-eval optimum 0.1928\n', bestJ);
fprintf('  saved: %s\n', fullfile(RUNDIR,'longfit_hist.mat'));
fprintf('  HOLD: slope ~0.413 | level_c in [3.5 6.5] | amp >= 40 | osc <= -40\n');
fprintf('  WIN:  maperr below 80 (project best ever 79.17; pre-today record 90.97)\n');
disp('LONGFIT_DONE');
