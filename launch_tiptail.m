% Fit the 2-CHAMBER for tip-tail contrast / latency exponent d.
% Native 2-chamber already gives d=0.445 (closest any config has come to the
% 0.39-0.41 target). Nothing has ever been fitted FOR tip-tail contrast, so this
% is a new objective, not a re-run. Includes k3e/r3e/k4e (17-19) -- the PLACE
% DEPENDENCE of the active terms -- which the default fitidx omits and which is
% exactly the lever for a tip that collapses at high CF.
pa=modpar26(2);
o.warm    = pa;
o.tiptail = 1;  o.dtarget = 0.40;  o.wd = 1;
o.fitidx  = [1 2 3 7 8 9 10 11 12 13 17 18 19 20];  % k1o r1o m1o k3o r3o k4o aco
                                                    % k1e r1e m1e k3e r3e k4e ace
o.maptol  = 185;  o.wmap = 0.02;
o.wcrit   = 0.05;
o.maxfe   = 300;                 % affordable: ~30 s/eval
o.out     = 'parfit26_tiptail2.mat';
R = parfit26(2,o);
m0=tiptail_metric(modpar26(2)); mf=tiptail_metric(R.pa);
fprintf('\n=== TIP-TAIL FIT (2-chamber) ===\n');
fprintf('  d       : %.3f -> %.3f   (target 0.39-0.41)\n', m0.d, mf.d);
fprintf('  hiCF con: %.1f -> %.1f dB\n', m0.chi, mf.chi);
fprintf('  nvalid  : %d -> %d      R2 %.2f -> %.2f\n', m0.nvalid, mf.nvalid, m0.r2, mf.r2);
fprintf('  maperr  : %.1f (tol 185)   osc %+.1f\n', R.Rf.maperr, R.S.maxRe_osc);
fprintf('  BF seq  : %s\n', num2str(mf.BF,'%7.2f'));
fprintf('  contrast: %s\n', num2str(mf.contrast,'%7.1f'));
disp('TIPTAIL_FIT_DONE');
