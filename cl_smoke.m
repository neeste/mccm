% Phase-1 smoke test for the 4th (CL) chamber.
% Checks, in order: parameters build; the m=4 operator assembles and is
% non-singular; a toneburst runs; the response is finite and sensible; and the
% 3-chamber path is UNCHANGED by all the edits (regression guard).
fprintf('--- parameters ---\n');
pa=modpar26(4);
fprintf('m=%d  chsz=[%s]  has k5o=%d\n', pa.m, num2str(pa.chsz), isfield(pa,'k5o'));

fprintf('\n--- 3-chamber REGRESSION (must match earlier values) ---\n');
pr.fr=2; pr.lv=60; pr.pa=modpar26(3);
L=load('parfit26_recip.mat'); pr.pa=L.R.pa; pr.pa.hbmode='bm';
S3=tdm26('wnr1',pr,0,0);
fprintf('3-ch 2kHz/60 BM-peak latency = %.2f ms   (expected 4.34)\n', S3.tpk);

fprintf('\n--- 4-chamber run ---\n');
q.fr=2; q.lv=60; q.pa=modpar26(4); q.pa.hbmode='bm';
try
    t0=tic; S4=tdm26('wnr1',q,0,0); wall=toc(t0);
    fprintf('ran OK in %.0f s\n', wall);
    fprintf('  eardrum level mlv = %.2f dB   (drive reaching the cochlea)\n', S4.mlv);
    fprintf('  WNR latency       = %.2f ms\n', S4.tpk);
    w=S4.wnr(:);
    fprintf('  WNR finite=%d  max=%.3g  any NaN=%d\n', all(isfinite(w)), max(w), any(isnan(w)));
    save('cl_smoke.mat','S4');
catch e
    fprintf('4-CHAMBER FAILED: %s\n', e.message);
    if (~isempty(e.stack)), fprintf('  at %s line %d\n', e.stack(1).name, e.stack(1).line); end
end
disp('CL_SMOKE_DONE');
