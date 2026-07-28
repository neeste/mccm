% SMOKE TEST: does parfit26 actually work for nch=4, and does J MOVE?
% parfit26 has never been run with 4 chambers. Before spending a long fit,
% confirm (a) it does not throw, (b) J is not pinned at 1e6 (which is what a
% silently-failing fdm26 would produce -- the exact failure this whole fdm26
% extension was about), and (c) J CHANGES between evaluations.
%
% Pure MAP fit: all ABR-dependent weights zero, so parfit26 skips abr_metric and
% each eval is just fdm26 (seconds). cheapstab=1 skips the per-eval coupeig,
% which for m=4 (N=8406, dense eig) would be ~4 min/eval.

opts = struct();
opts.tiptail   = 0;          % no ABR terms at all -> pure tuning/MAP objective
opts.wslope    = 0; opts.wlevel = 0; opts.wanchor = 0;
opts.wshoulder = 0; opts.wlcb   = 0; opts.surface = 0;
opts.wmap      = 0.01;       % the only active term
opts.maptol    = 100;        % ambitious floor (m=2 native is 104.6), so the
                             % maperr term stays active for the whole fit
opts.wcrit     = 0.005;
opts.cheapstab = 1;          % skip per-eval coupeig; verify at the end
opts.warm      = 'none';     % no warm start (refit_c3_map.mat is 3-chamber)
opts.hbmode    = 'bm';
% CF-map parameters: k1o m1o aco / k1e m1e ace / k1q m1q / chsz(1:4).
% The tonotopic RANGE is set mainly by the exponential slopes k1e,m1e.
opts.fitidx    = [1 3 10 11 13 20 21 23 31 32 33 34];
opts.maxfe     = 6;          % smoke only
opts.out       = 'smoke_c4map.mat';

t0=tic; R = parfit26(4, opts); w=toc(t0);
fprintf('\n=== SMOKE RESULT (%.0f s for %d evals) ===\n', w, opts.maxfe);
if (isstruct(R))
    fn=fieldnames(R);
    fprintf('  returned fields: %s\n', strjoin(fn(1:min(8,numel(fn)))', ', '));
    if (isfield(R,'Rf') && isfield(R.Rf,'maperr'))
        fprintf('  final maperr : %.2f   (native was 1022.96)\n', R.Rf.maperr);
    end
    if (isfield(R,'S'))
        if (isfield(R.S,'maxRe')),     fprintf('  maxRe(all)   : %+.1f\n', R.S.maxRe); end
        if (isfield(R.S,'maxRe_osc')), fprintf('  maxRe_osc    : %+.1f\n', R.S.maxRe_osc); end
    end
end
fprintf(['\nPASS if it ran without throwing AND maperr is finite and not stuck at\n' ...
         'the native 1022.96 (J pinned at 1e6 would leave it unchanged).\n']);
disp('SMOKE_C4MAP_DONE');
