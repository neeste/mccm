% 4-CHAMBER MAP REFIT, HEAVY m2o PINNED, WITH A tdm26 VALIDITY GUARD.
%
% The first heavy-m2 refit (parfit26_c4map_m2) was a pure-maperr fit (tiptail=0,
% cheapstab=1) with NO tdm26-side constraint. It reached maperr 464 in fdm26 but
% the tdm26 click DIVERGED (12852/14336 non-finite): the fit crushed SV to 0.37
% (native 1.0), a stiff region fdm26's steady-state solve and coupeig's
% eigenvalues both bless but the time march cannot survive.
%
% FIX (no code change): tiptail=1 makes jointobj call tiptail_metric every eval,
% which time-marches the tdm26 click and returns J=1e6 the instant it diverges
% or loses its CF map (the tt.ok gate, parfit26.m:146). wd=0 so it does NOT chase
% a d target -- it is purely the validity guard -- while wmap keeps maperr the
% objective and the built-in chi/nvalid terms push toward a sharp, intact click
% map. cheapstab=1 + tiptail=1 is the CORRECT pairing (parfit26.m:157-160).
% Cost ~2x/eval (fdm26 + tdm26 click), ~40 s, so ~3 h for 300 evals.

b2 = modpar26(4).m2o;
opts = struct();
opts.pin       = struct('m2o', b2*32);
opts.tiptail   = 1;        % ON = tdm26-click validity guard via tt.ok
opts.wd        = 0;        % but do NOT chase d -- guard only, maperr is the goal
opts.wmap      = 0.01;
opts.maptol    = 100;
opts.wcrit     = 0.005;
opts.cheapstab = 1;        % correct with tiptail=1 (tt screen is the stability proxy)
opts.skipabr   = 1;        % m=4 abr_metric ~38 min/end; not needed for a map fit
opts.warm      = 'none';
opts.hbmode    = 'bm';
opts.fitidx    = [1 3 10 11 13 20 21 23 31 32 33 34];  % k1o m1o aco / k1e m1e ace / k1q m1q / chsz1-4
opts.maxfe     = 300;
opts.out       = 'parfit26_c4map_m2b.mat';

fprintf('pinned m2o=%.4g (x32);  tiptail guard ON, wd=0, wmap=%.3g\n', b2*32, opts.wmap);
t0=tic; R = parfit26(4, opts); w=toc(t0);
fprintf('\n=== C4 MAP-m2b FIT RESULT (%.0f s) ===\n', w);
if (isstruct(R))
    if (isfield(R,'Rf')&&isfield(R.Rf,'maperr')), fprintf('  final maperr : %.2f  (broken m2-fit was 464 but DIVERGED; c4map 547.8)\n', R.Rf.maperr); end
    if (isfield(R,'pa')),  fprintf('  pin held? R.pa.m2o = %.4g (want %.4g);  chsz=[%s]\n', R.pa.m2o, b2*32, num2str(R.pa.chsz,'%.2f ')); end
    if (isfield(R,'S')&&isfield(R.S,'maxRe')), fprintf('  maxRe(all) %+.1f  osc %+.1f\n', R.S.maxRe, R.S.maxRe_osc); end
end
disp('C4MAP_M2B_DONE');
