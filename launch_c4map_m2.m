% 4-CHAMBER MAP REFIT WITH HEAVY m2o PINNED.
%
% The first map fit (parfit26_c4map) could not decompress the CF map because its
% fitidx had chsz but NOT the mass couplings -- and the mass sweep showed m2o is
% the lever that moves RANGE (native 3.43 -> 4.68 oct at m2o x32, fold removed,
% amplifier intact +2.46 dB, stable). So pin m2o heavy and refit the impedance/
% area profiles around it, to (a) close more of the 1.2-oct gap to 5.90 and
% (b) recover some of the contrast lost to the heavy mass (13.5 -> 9.4 dB).
%
% Pin uses the staged opts.pin machinery: m2o is absent from parnames, so
% setpar_l never touches it and it persists on base + warm start. Pure MAP fit
% (skipabr=1, cheapstab=1). coupeig still runs in the FINAL reconstruction, so
% end-stability is verified even without a per-eval guard.

b2 = modpar26(4).m2o;
opts = struct();
opts.pin       = struct('m2o', b2*32);   % heavy m2o (the best-range config)
opts.tiptail   = 0;
opts.wslope=0; opts.wlevel=0; opts.wanchor=0; opts.wshoulder=0; opts.wlcb=0; opts.surface=0;
opts.wmap      = 0.01;
opts.maptol    = 100;
opts.wcrit     = 0.005;
opts.cheapstab = 1;
opts.skipabr   = 1;
opts.warm      = 'none';                  % start from native profiles
opts.hbmode    = 'bm';
opts.fitidx    = [1 3 10 11 13 20 21 23 31 32 33 34];  % k1o m1o aco / k1e m1e ace / k1q m1q / chsz1-4
opts.maxfe     = 300;
opts.out       = 'parfit26_c4map_m2.mat';

fprintf('pinned m2o = %.4g (x32 of native %.4g)\n', b2*32, b2);
t0=tic; R = parfit26(4, opts); w=toc(t0);
fprintf('\n=== C4 MAP-m2 FIT RESULT (%.0f s) ===\n', w);
if (isstruct(R))
    if (isfield(R,'Rf')&&isfield(R.Rf,'maperr')), fprintf('  final maperr : %.2f  (native-m2 start below; c4map was 547.8)\n', R.Rf.maperr); end
    if (isfield(R,'pa')),  fprintf('  pinned m2o held? R.pa.m2o = %.4g (want %.4g)\n', R.pa.m2o, b2*32); end
    if (isfield(R,'S')&&isfield(R.S,'maxRe')), fprintf('  maxRe(all)   : %+.1f   maxRe_osc %+.1f\n', R.S.maxRe, R.S.maxRe_osc); end
end
disp('C4MAP_M2_DONE');
