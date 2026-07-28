% TBOAE run: 3-chamber, roughness above floor, OAE as secondary TBABR output.
pa=modpar26(3); pa.oae=1; pa.rough_amp=3e-2;
S=tdm26('tbabr',pa,0,0);
r=S.oae./S.lat;                         % OAE latency / WNR(forward) latency
fprintf('\n=== TBOAE smoke (3-chamber, rough_amp=1e-4) ===\n');
fprintf(' f(kHz) lvl  WNR(ms) OAE(ms) ratio  emag(dB)\n');
for j=1:4
  for k=1:4
    fprintf('%6.2f %3.0f  %6.2f  %6.2f  %5.2f  %6.1f\n', ...
        S.f(j), S.slv(k), S.lat(j,k), S.oae(j,k), r(j,k), S.oam(j,k));
  end
end
fprintf('\nOAE/WNR ratio, mean over levels (2013: ~2x above 1.5kHz, ~1.3x below):\n');
for j=1:4, fprintf('  %.2f kHz : %.2f\n', S.f(j), mean(r(j,:),'omitnan')); end
fprintf('\nfinite OAE cells: %d/16   median emag: %.1f dB\n', ...
        sum(isfinite(S.oae(:))), median(S.oam(isfinite(S.oam))));
save('tboae_smoke.mat','S');
disp('TBOAE_SMOKE_DONE');
