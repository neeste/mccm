% How expensive is fdm26 for m=3? Decides whether an amplifier term in the
% objective (needs active + passive runs) is affordable per evaluation.
pa=modpar26(3);
t0=tic; R1=fdm26(struct('pa',pa)); t1=toc(t0);
pa0=pa; pa0.gam=0;
t0=tic; R0=fdm26(struct('pa',pa0)); t2=toc(t0);
fprintf('\n fdm26 m=3 active : %.1f s   maperr %.1f  tipgain %.1f dB\n', t1, R1.maperr, R1.tipgain);
fprintf(' fdm26 m=3 passive: %.1f s   maperr %.1f  tipgain %.1f dB\n', t2, R0.maperr, R0.tipgain);
fprintf(' fd amp proxy (active tipgain - passive tipgain) = %+.2f dB\n', R1.tipgain-R0.tipgain);
fprintf(' [score26 click-based amp for this config was +81.15 dB]\n');
fprintf(' => 2 fdm26 calls per eval costs %.1f s; 300 evals = %.1f min\n', t1+t2, 300*(t1+t2)/60);
disp('TIME_FDM_DONE');
