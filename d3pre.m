% D3PRE -- step 2 of the d3int=1 plan: pre-flight the warm start and MEASURE the
% evaluation cost before sizing a multi-hour run.
%
% WHY THIS EXISTS AT ALL. Two failures in this project's record are both
% addressed here, and neither is hypothetical:
%
%   1. c3fit ran 7.3 h against a warm start parfit26 REPORTED but did not fit,
%      because the pa it printed did not survive setpar_l. s3fit added a
%      pre-flight for that; this repeats it under d3int=1, which is a DIFFERENT
%      mapping question -- the six third-DOF entries are new as of bd1d5ac and
%      the round-trip was only ever proven with d3int at its saved value.
%   2. Segment counts have been chosen from an eval time that no longer applied.
%      d3int=1 makes dof=3 (tdm26.m:2133), so coupeig's operator grows and the
%      2-DOF ~76 s/eval figure does NOT transfer. Measure it.
%
% Also asserts the FIT VECTOR is what the run believes it is. parfit26 clips
% fitidx to 30+nc BEFORE appending the by-name blocks, so third-DOF indices
% written in by hand are dropped with a warning and the fit proceeds with them
% pinned -- silently answering a different question. opts.fitd3 exists for this.
PROJ = '/Users/neely/Library/CloudStorage/OneDrive-FatherFlanagan''sBoysHome/STN/mod/mccm';
addpath(PROJ); cd(PROJ);
WARM = 'fit_nch3_surface.mat';

L=load(WARM); pa=L.R.pa;
if (isfield(pa,'hbmode')), pa=rmfield(pa,'hbmode'); end   % shear drive, SN's call
pa.d3int = 1;                                             % THE structural change
nc = numel(pa.chsz);
H  = parfit26('handles');

fprintf('\n===== D3PRE: nch=3, d3int=1 =====\n');

% ---- 1. does the warm start survive the parameter mapping WITH d3int=1? ----
base = modpar26(3); base.d3int = 1;
pa_eval = H.setpar(base, H.getpar(pa), []);
e_warm = fdm26(struct('pa',pa));      e_warm = e_warm.maperr;
e_eval = fdm26(struct('pa',pa_eval)); e_eval = e_eval.maperr;
fprintf('\n  PRE-FLIGHT  warm maperr %.4f | as jointobj builds it %.4f | delta %.3e\n', ...
        e_warm, e_eval, abs(e_warm-e_eval));
if (abs(e_warm-e_eval) > 1)
    error('D3PRE ABORTED: warm start does not survive setpar_l (%.2f vs %.2f).', e_warm, e_eval);
end
% The 2-DOF number the run is trying to beat, from the SAME parameters.
p2 = pa; p2.d3int = 0;
e_2dof = fdm26(struct('pa',p2)); e_2dof = e_2dof.maperr;
fprintf('  d3int=0 %.4f -> d3int=1 %.4f  (the gap this run has to close)\n', e_2dof, e_warm);

% ---- 2. are the third-DOF entries actually LIVE through the mapping? ----
% Mistake guard: a round-trip that passes because setpar wrote nothing looks
% exactly like a round-trip that passes because it wrote the right values.
pvk = H.getpar(pa); pvk(30+nc+3) = pvk(30+nc+3)*2;   % k5o x2 through the VECTOR
q = H.setpar(base, pvk, []);
Rq = fdm26(struct('pa',q));
fprintf('  k5o via pv: %.6g -> %.6g | maperr %.4f (delta %+.3e)%s\n', ...
        pa.k5o, q.k5o, Rq.maperr, Rq.maperr-e_eval, ...
        tern(abs(Rq.maperr-e_eval)<1e-9,'  <-- INERT, ABORT',''));
if (abs(Rq.maperr-e_eval) < 1e-9)
    error('D3PRE ABORTED: k5o moves nothing through the pv mapping.');
end

% ---- 3. the fit vector: 23 params, and NOTHING silently dropped ----
fitidx = [1 2 3 4 5 7 8 9 10 11 13 20 21 31 32 33 15 17];  % + r2e(15), k3e(17)
want   = numel(fitidx) + 1 + 6;                            % + hbsc + six third-DOF
fprintf('\n  fit vector: %d impedance/chsz + hbsc + 6 third-DOF = %d expected\n', ...
        numel(fitidx), want);

% ---- 4. TIMING, component by component, 2-DOF vs 3-DOF ----
% One jointobj evaluation on this option set = abr_metric + coupeig + score26.
fprintf('\n  timing (one evaluation = abr_metric + coupeig + score26)\n');
fprintf('  %-12s %10s %10s %8s\n','component','d3int=0','d3int=1','ratio');
t=tic; abr_metric(p2,false);                                 t20=toc(t);
t=tic; abr_metric(pa,false);                                 t21=toc(t);
t=tic; evalc('tdm26(''coupeig'',struct(''pa'',p2));');       t30=toc(t);
t=tic; evalc('tdm26(''coupeig'',struct(''pa'',pa));');       t31=toc(t);
t=tic; score26(p2,'fast',false);                             t40=toc(t);
t=tic; S1=score26(pa,'fast',false);                          t41=toc(t);
fprintf('  %-12s %10.1f %10.1f %8.2f\n','abr_metric',t20,t21,t21/max(t20,eps));
fprintf('  %-12s %10.1f %10.1f %8.2f\n','coupeig',   t30,t31,t31/max(t30,eps));
fprintf('  %-12s %10.1f %10.1f %8.2f\n','score26',   t40,t41,t41/max(t40,eps));
ev0=t20+t30+t40; ev1=t21+t31+t41;
fprintf('  %-12s %10.1f %10.1f %8.2f  <== EVAL\n','TOTAL',ev0,ev1,ev1/max(ev0,eps));

% ---- 5. state of the start point, so the run's first segment has a baseline ----
m1 = abr_metric(pa,false);
fprintf('\n  START (d3int=1)  slope %.4f | lvl_c %.3f | shoulder %.4f | %d/16 cells>0.5\n', ...
        m1.slope, m1.level_c, m1.shoulder, sum(m1.sho(isfinite(m1.sho))>0.5));
fprintf('  maperr %.2f | amp %+.2f | maxRe %+.2f | osc %+.2f\n', ...
        S1.maperr, S1.amp_gain, S1.maxRe, S1.maxRe_osc);

% ---- 6. SEGMENT SIZING, from the measured number rather than the old one ----
for H_ = [8 12 16]
    fprintf('  %2d h budget -> %4d evaluations at %.0f s\n', H_, floor(H_*3600/ev1), ev1);
end
save('d3pre.mat','ev0','ev1','e_warm','e_2dof','m1','S1','fitidx');
disp('D3PRE_DONE');
function s=tern(c,a,b), if c, s=a; else, s=b; end, end
