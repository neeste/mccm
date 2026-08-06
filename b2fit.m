% B2FIT -- push the BEST nch=3 state further, targeting its one remaining defect.
%
% ======================== WHY THIS AND NOT AN m=4 FIT ====================
%
% m=4 is the more interesting unknown (amp +54.06 dB, maperr 329.32, lvl_c in
% band) but it is not affordable here. Three calls at m=4 exceeded 10 minutes
% while sharing CPU with rlfit, against 135 s for the same three at m=3. Even
% discounting contention that is >= 300 s/eval, giving ~100 evaluations in the
% window for a 19-parameter vector -- below the 20-point simplex plus useful
% search. m=4 needs its own dedicated window and a measured eval time first.
%
% ==================== WHAT THIS RUN IS ACTUALLY AFTER ====================
%
% d3fit arm B is the best three-chamber state on record:
%
%   slope 0.4899 (target 0.413) | lvl_c 5.672 IN BAND | maperr 147.92 IN
%   THRESHOLD | amp +43.34 dB | osc -170.3
%
% and it has exactly ONE defect left: the second peak. Its mean shoulder is
% 0.2068, which is 12 cells at ~0.000 and FOUR cells at 0.663-0.992, all at 60
% and 80 dB SPL, one of them a coin flip at 0.992 (dsweep.m).
%
% Those four cells are why every latency this project quotes still carries a
% caveat. The detector can take either peak there, so the reported latency in
% those cells is not a measurement of the model. Driving them to zero would make
% the nch=3 latency surface quotable WITHOUT the caveat for the first time, and
% that is worth more to MOH2027 than another decimal place on the slope.
%
% IT IS STILL MOVING, which is what makes this worth budget rather than hope.
% Arm B's two segments went J 1.0615 -> 1.0566 (0.46%) -> 0.9633 (8.8%), and
% nearly all of the second segment was the shoulder coming down 0.2603 ->
% 0.2068. A converged run does not do that on its last segment.
%
% WEIGHTS UNCHANGED from d3fit arm B, deliberately. Raising wshoulder would
% attack the four cells harder but would make this a different objective, so the
% result could not be appended to arm B's curve and the comparison with the
% three-degree-of-freedom arm would break. Same objective, more evaluations.
% If the shoulder stalls here, THEN a reweighted run is the next experiment and
% it will have this curve to be judged against.
% =======================================================================
RUNDIR = '/Users/neely/mccm_runs';
if (~exist(RUNDIR,'dir')), mkdir(RUNDIR); end
diary(fullfile(RUNDIR,'b2fit.log')); diary on;
cln = onCleanup(@() diary('off'));

NSEG = 2; FE = 150;
WSHOULDER = 2.0; WSURF = 1.0; HBMODE = ''; WGAIN = 0.01;
FITIDX = [1 2 3 4 5 7 8 9 10 11 13 20 21 31 32 33 15 17];
PIN = struct('d3int',0);            % the CONTROL configuration: 2-DOF

WARM = fullfile(RUNDIR,'d3fit_B_seg2.mat');
if (~exist(WARM,'file')), error('b2fit ABORTED: %s not found.', WARM); end
L=load(WARM); pa=L.R.pa;
if (isfield(pa,'hbmode')), pa=rmfield(pa,'hbmode'); end
pa.d3int = 0;
if (isfield(pa,'rlsplit')), pa.rlsplit = 0; end   % never on for this arm

if (isempty(HBMODE) && isfield(pa,'hbmode') && ~isempty(pa.hbmode))
    error('b2fit ABORTED: HBMODE is shear but warm start carries hbmode="%s".', pa.hbmode);
end
if (isfield(pa,'latsoft') && isinf(pa.latsoft))
    error('b2fit ABORTED: warm start sets latsoft=Inf (legacy gameable detector).');
end

fprintf('\n===== B2FIT: continue the best nch=3 state, %d x %d evals =====\n', NSEG, FE);
fprintf('  warm: d3fit arm B seg2, J 0.9633\n');
fprintf('  TARGET: the 4 cells at 60-80 dB with shoulder 0.663-0.992\n');
fprintf('  weights UNCHANGED from arm B so this extends its curve\n');

H=parfit26('handles'); b0=modpar26(3); b0.d3int=0;
pe=H.setpar(b0,H.getpar(pa),[]);
ew=fdm26(struct('pa',pa)); ew=ew.maperr;
ee=fdm26(struct('pa',pe)); ee=ee.maperr;
fprintf('  PRE-FLIGHT  warm %.4f | rebuilt %.4f | delta %.3e\n', ew, ee, abs(ew-ee));
if (abs(ew-ee) > 1), error('b2fit ABORTED: warm start does not survive setpar_l.'); end

m0 = abr_metric(pa,false);
fprintf('  START shoulder %.4f | %d/16 cells>0.5 | max %.3f\n', ...
        m0.shoulder, sum(m0.sho(isfinite(m0.sho))>0.5), max(m0.sho(isfinite(m0.sho))));

hist=struct('seg',{},'J',{},'slope',{},'lvlc',{},'sho',{},'ncell',{},'maperr',{},'amp',{},'osc',{},'hrs',{});
T0=tic;
for s=1:NSEG
    t0=tic; out=fullfile(RUNDIR,sprintf('b2fit_seg%d.mat',s));
    fprintf('----- segment %d/%d -----\n', s, NSEG);
    try
        R = parfit26(3, struct('maxfe',FE, 'surface',1, 'wsurf',WSURF, ...
                               'wshoulder',WSHOULDER, 'hbmode',HBMODE, ...
                               'wgain',WGAIN, 'fithbsc',true, 'fitidx',FITIDX, ...
                               'pin',PIN, 'warm',pa, 'verbterm',true, 'out',out));
    catch e
        fprintf('  SEGMENT %d FAILED: %s\n', s, e.message(1:min(200,end))); break;
    end
    pa = R.pa; S = score26(pa,'fast',false); m = abr_metric(pa,false);
    nc5 = sum(m.sho(isfinite(m.sho))>0.5);
    hist(end+1)=struct('seg',s,'J',R.J,'slope',m.slope,'lvlc',m.level_c, ...
        'sho',m.shoulder,'ncell',nc5,'maperr',S.maperr,'amp',S.amp_gain, ...
        'osc',S.maxRe_osc,'hrs',toc(t0)/3600); %#ok<SAGROW>
    fprintf('\n  B2 SEG %d | J %.4f | SHOULDER %.4f (%d/16 cells>0.5) | slope %.4f\n', ...
            s, R.J, m.shoulder, nc5, m.slope);
    fprintf('           | lvl_c %.3f | maperr %.2f | amp %+.2f | osc %.1f | %.1f h\n', ...
            m.level_c, S.maperr, S.amp_gain, S.maxRe_osc, hist(end).hrs);
    save(fullfile(RUNDIR,'b2fit_hist.mat'),'hist');
end

fprintf('\n===== B2FIT DONE, %.1f h =====\n', toc(T0)/3600);
fprintf('  seg |      J | shldr  | >0.5 | slope  | lvl_c | maperr |   amp  | hrs\n');
for i=1:numel(hist)
    h=hist(i);
    fprintf('  %3d | %6.4f | %6.4f | %4d | %6.4f | %5.3f | %6.2f | %+6.2f | %4.1f\n', ...
            h.seg, h.J, h.sho, h.ncell, h.slope, h.lvlc, h.maperr, h.amp, h.hrs);
end
fprintf('\n  READ THE >0.5 COLUMN. It was 4 of 16 at the start. Reaching 0 is what\n');
fprintf('  would let the nch=3 latency surface be quoted without a caveat; the mean\n');
fprintf('  shoulder falling while that count holds at 4 is NOT the same result.\n');
disp('B2FIT_DONE');
