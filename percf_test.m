function percf_test
% Per-CF roughness window (DIAGNOSTIC probe): roughness confined to +-0.75 oct of
% each stimulus CF, so the CF coherent reflection is the emission. Does the
% measured OAE latency now track 2xWNR (the 2013 high-frequency signature)?
fprintf('\n=== per-CF roughness window : does OAE track 2xWNR? (60 dB) ===\n');
fprintf(' f(kHz)  WNR    OAE   2xWNR   OAE/WNR  emag(dB)\n');
for fr=[0.5 1 2 4]
  pr.fr=fr; pr.lv=60; pr.pa=modpar26(3);
  pr.pa.oae=1; pr.pa.rough_amp=3e-2; pr.pa.rough_percf=1;
  S=tdm26('wnr1',pr,0,0);
  fprintf('%6.2f  %5.2f  %5.2f  %5.2f    %4.2f    %6.1f\n', ...
          fr, S.tpk, S.oae, 2*S.tpk, S.oae/S.tpk, S.oam);
end
disp('PERCF_TEST_DONE');
