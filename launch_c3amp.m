% AMPLIFIER-PRESERVING 3-CHAMBER MAP FIT.
%
% The native 3-chamber has amp +81.15 dB -- the strongest amplifier measured
% anywhere in this project (m=1 champion +50.68, m=2 native +39.11) -- with
% maperr 499.3. The fitted champion refit_c3_map has maperr 164 but amp only
% +40: THE FIT TRADED AWAY 41 dB OF AMPLIFIER, invisibly, because no objective
% ever contained an amplifier term and no scorecard measured one before score26.
%
% QUESTION: how much maperr can be recovered WITHOUT giving up the gain?
%
% METHOD -- preserve the amplifier BY CONSTRUCTION rather than by measuring it.
% Measuring active-vs-passive gain per evaluation needs two extra model runs
% (~4 h for 300 evals), and fdm26's cheap tipgain proxy is NOT a substitute: it
% reads +23.70 dB where the validated click measurement reads +81.15 (different
% observable -- tip-to-tail on Dh vs active-vs-passive on d1). So instead:
%   FREEZE the ACTIVE impedances k3o(7), r3o(8), k4o(9) at their native values
%   and fit only the passive/tuning parameters.
% k_act = gh*k3 - gam*k4 and r_act = gh*r3 - gam*r4, so holding k3/r3/k4 fixed
% (with gam=1 and gpo/gpe untouched) keeps the active force intact while the
% tuning/CF-map parameters move freely. Zero extra cost, and directly testable
% afterwards with score26.
%
% fitidx: [1 2 3 4 5 10 11 13 20 21 31 32 33]
%         k1o r1o m1o k2o r2o aco | k1e m1e ace k1q | chsz1-3
%         = the DEFAULT fitidx MINUS 7,8,9 (k3o r3o k4o) -- the active terms.
%
% GUARD: tiptail=1 with wd=0 gives a per-eval tdm26 validity check (tt.ok
% rejects a diverged or CF-map-broken model). This works for m=3 (nvalid=5) --
% it is only m=4 where tiptail_metric's 0.6-16 kHz window fails. Its built-in
% chi term also mildly rewards keeping the tip. statol is now 100 (calibrated).
%
% PRE-FLIGHT FIRST -- the hard rule after the 7-hour loss: evaluate the guard at
% the START POINT and confirm it passes before spending the run.

fprintf('=== PRE-FLIGHT: does the tdm26 validity guard pass at the start point? ===\n');
pa0 = modpar26(3);
t = tiptail_metric(pa0, false);
fprintf('  native m=3: tt.ok=%d  nvalid=%d  d=%.3f  R2=%.2f  %s\n', ...
        t.ok, t.nvalid, t.d, t.r2, t.msg);
if (~(t.ok || t.nvalid >= 4))
    fprintf('  *** ABORT: the guard rejects the START POINT -- every eval would be\n');
    fprintf('  J=1e6 and the fit would be flat (the 7-hour failure mode). ***\n');
    disp('C3AMP_ABORTED'); return
end
fprintf('  guard passes -- proceeding.\n\n');

opts = struct();
opts.warm      = 'none';        % start from NATIVE m=3 (which HAS the +81 dB),
                                % not refit_c3_map (which already lost it)
opts.fitidx    = [1 2 3 4 5 10 11 13 20 21 31 32 33];
opts.tiptail   = 1;             % validity guard only
opts.wd        = 0;             % do NOT chase d (it is noise-dominated anyway)
opts.wmap      = 0.01;
opts.maptol    = 100;
opts.wcrit     = 0.005;
opts.statol    = 100;           % calibrated: healthy models reach +23.6
opts.cheapstab = 1;             % per-eval coupeig too slow for m=3 (~86 s)
opts.skipabr   = 1;
opts.hbmode    = 'bm';
opts.maxfe     = 300;
opts.out       = 'parfit26_c3amp.mat';

t0=tic; R = parfit26(3, opts); w=toc(t0);
fprintf('\n=== C3AMP FIT (%.0f s) ===\n', w);
if (isstruct(R) && isfield(R,'Rf')), fprintf('  maperr %.1f  (native 499.3, refit_c3_map 164.0)\n', R.Rf.maperr); end

% ---- the decisive check: did the amplifier survive? ------------------------
fprintf('\n=== VERIFY with score26 (did the +81 dB survive?) ===\n');
S = score26(R.pa, 'fast');
fprintf('\n  amp %+.2f dB (native +81.15, refit_c3_map +40.13)  maperr %.1f\n', ...
        S.amp_gain, S.maperr);
save('parfit26_c3amp_score.mat','S');
disp('C3AMP_DONE');
