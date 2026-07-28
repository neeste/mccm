% SHAPE FIT OF THE 1-CHAMBER CHAMPION.
%
% WHY SHAPE, NOT d. The scalar d is NOISE-DOMINATED: per-level slopes scatter
% 0.36-0.38 across 20-80 dB, WIDER than the 0.6->0.4 distance to target, so
% every fit that chased d to three decimals was partly fitting noise. The SHAPE
% residual -- rms log-residual of all 16 latency cells against
% tau = b*c^(-i)*f^(-d) -- is 0.132 for this champion: a 13% structural error
% that is large, well-determined, and never targeted.
%
% m=1 is the right subject: it is the best model of the TBABR+TBOAE data on
% every other axis (maperr 74.6, b=12.62 IN BAND, level 1.61 %/dB vs target
% 1.62, OAE ratio 2.25, clean 6.37-oct map, +50.7 dB amplifier). Its ONLY miss
% is the latency-surface shape. Improving that improves the best model rather
% than re-running another chamber-count experiment.
%
% OBJECTIVE: wshape * abr_surface_obj J  (= shape-rms + hinge losses on b,c,d
% that are ZERO inside the published uncertainty, so the fit is not driven to
% spurious precision) + wmap * max(0, maperr - maptol) to PROTECT the tuning
% that already works. maptol=80 keeps maperr near the champion's 74.6.
%
% NOTE fitidx EXCLUDES index 33: for m=1, chsz has only TWO elements, so
% pv has length 30+2 and index 33 does not exist (the default fitidx assumes 3).

fprintf('=== PRE-FLIGHT: evaluate the objective chain at the START point ===\n');
L = load('refit_c1broad.mat'); pa0 = L.R.pa; pa0.hbmode = 'bm';
t0 = tic; m0 = abr_metric(pa0, false); tm = toc(t0);
fprintf('  abr_metric: ok=%d  n_sub=%d  (%.0f s/eval)\n', m0.ok, m0.n_sub, tm);
if (~m0.ok)
    fprintf('  *** ABORT: abr_metric fails at the start point. ***\n'); disp('C1SHAPE_ABORTED'); return
end
[Js, Ds] = abr_surface_obj(m0);
fprintf('  abr_surface_obj: J=%.4f  shape-rms=%.4f  b=%.2f c=%.2f d=%.3f  n=%d\n', ...
        Js, Ds.resid, Ds.b, Ds.c, Ds.d, Ds.n);
fprintf('  hinges d/c/b = %.3f / %.3f / %.3f\n', Ds.hd, Ds.hc, Ds.hb);
if (~isfinite(Js) || Js >= 1e5)
    fprintf('  *** ABORT: objective is non-finite or at the sentinel. ***\n'); disp('C1SHAPE_ABORTED'); return
end
fprintf('  objective is LIVE (J=%.4f). Estimated %.0f min for 150 evals.\n\n', Js, 150*tm/60);

opts = struct();
opts.warm      = pa0;           % start FROM the champion, so any gain is banked on top
opts.fitidx    = [1 2 3 4 5 7 8 9 10 11 13 20 21 31 32];   % NOTE: no 33 (m=1 has 2 chsz)
opts.wshape    = 1;             % THE objective
opts.wslope=0; opts.wlevel=0; opts.wanchor=0; opts.wshoulder=0; opts.wlcb=0;
opts.surface   = 0;             % NOT the surf_term band objective -- the shape one
opts.tiptail   = 0;
opts.wmap      = 0.02;          % protect the tuning that already works
opts.maptol    = 80;            % champion is 74.6
opts.wcrit     = 0.005;
opts.statol    = 100;           % calibrated
opts.cheapstab = 0;             % m=1 coupeig is cheap enough to keep per-eval
opts.skipabr   = 0;             % abr_metric IS the objective here
opts.hbmode    = 'bm';
opts.maxfe     = 150;
opts.out       = 'parfit26_c1shape.mat';

tA=tic; R = parfit26(1, opts); w=toc(tA);
fprintf('\n=== C1SHAPE FIT (%.0f s) ===\n', w);

fprintf('\n=== VERIFY: full scorecard vs the champion ===\n');
S = score26(R.pa, 'full');
fprintf('\n  champion was: maperr 74.6 | shape 0.132 | b 12.62 | c 4.96 | d 0.609 | OAE 2.25 | amp +50.68\n');
save('parfit26_c1shape_score.mat','S');
disp('C1SHAPE_DONE');
