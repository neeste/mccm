% abr_run.m  -- Month-1 latency regression: the experiment sequence.
%
% Goal: recover the MOH24 Fig-3 result (WNR latency tracking ABR forward
% latency, slope d ~= 0.413 across frequency) that the MATLAB port lost at high
% frequencies, and map how far basal cochlear-amplifier gain can be pushed
% before the time-domain solver goes unstable.
%
% Run the blocks one at a time (Cmd+Enter per section), not all at once --
% each sweep is several tbabr runs and takes minutes.
%
% Reads:  tdm26.m, modpar26.m, abr_metric.m, abr_regress.m
% Writes: abr_baseline_*.mat, abr_sweep_*.mat, abr_gampro_*.mat  (+ figures)

%% 1) Baseline scorecard -- the number to beat.
%   Establishes today's error and, crucially, the model slope d. If d > 0.413
%   the high-frequency latencies are too short (the Fig-3 regression).
%   Expected wall time: nch=1 ~40 s. It prints a heartbeat before the run and a
%   scorecard after. If it runs for many minutes with no scorecard, it is NOT
%   computing (a healthy run is under a minute) -- interrupt with Ctrl-C and
%   check for a blocked figure window or a OneDrive file-lock on tbabr.txt.
R1 = abr_regress('baseline', 1);

%% 1b) 3-chamber baseline (~95 s).
%   With the wnr_latency onset detector the old edge-pinned artifacts (2 kHz/20,
%   4 kHz/20) are gone; any sub-threshold condition reports NaN and is excluded
%   (see m.n_sub). Expect slope ~0.59, level_c ~2.2 -- i.e. nch=3 does NOT fix
%   the latency slope (same problem as nch=1).
R3 = abr_regress('baseline', 3);

%% 2) Probe the active-coupling place profile: aco (basal magnitude).
%   cp.ac = aco*exp(ace*x); larger aco = more basal gain = more basal delay,
%   but also the first thing to destabilize. Watch for the OK->UNSTBL boundary.
Ra = abr_regress('sweep', 1, 'aco', [0.008 0.010 0.013 0.017 0.022 0.028]);

%% 3) Probe the active-coupling decay: ace (how basally weighted the gain is).
%   More negative ace concentrates gain at the base. Baseline is ~-0.40.
Re = abr_regress('sweep', 1, 'ace', [-0.30 -0.35 -0.40 -0.45 -0.50 -0.55]);

%% 4) Probe the named CA-gain profile pa.gampro directly: smooth basal boost.
%   gampro = 1 apically, (1+boost) at the stapes, ramped over basalfrac of L.
%   This is the cleanest "restore basal gain" knob and isolates the latency vs
%   instability trade-off you flagged. (In smoke testing, boost=2 already
%   diverged at basalfrac=0.35, so the stability ceiling is low -- these fine
%   steps bracket it; widen once you know where it sits.)
Rg = abr_regress('gampro', 1, [0 0.25 0.5 0.75 1.0 1.5], 0.35);

%% 5) (optional) Repeat the winning knob for the 3-chamber model.
Rg3 = abr_regress('gampro', 3, [0 0.5 1 2 4], 0.35);

% Interpreting the output:
%   - slope d dropping toward 0.413 as basal gain rises = latency recovering.
%   - the last OK value before UNSTBL = the stability ceiling on gain.
%   If d bottoms out ABOVE 0.413 while still stable, gain alone can't restore
%   Fig 3 and the fix belongs in the integrator (implicit/finer basal grid),
%   which is exactly the Month-2 decision this harness is meant to settle.
