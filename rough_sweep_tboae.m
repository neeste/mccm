% Diagnose the TBOAE emission vs roughness amplitude at one condition.
% Question 1: does emag scale linearly with rough_amp (real reflection) or
%             sit at a floor (numerical)?  +20 dB/decade == linear.
% Question 2: does the emission envelope PEAK near 2x the WNR forward latency
%             (physical round trip) or at ~record-center (~18 ms, noise)?
fr=2; lv=60;
ras=[1e-4 1e-3 1e-2 1e-1 3e-1];
fprintf('\n=== rough_amp sweep @ %g kHz, %g dB (3-chamber) ===\n',fr,lv);
fprintf('rough_amp   emag(dB)  centroid(ms)  envPeak(ms)  specGD(ms)   snr\n');
tpk=NaN; reclen=NaN;
for ra=ras
  pr.fr=fr; pr.lv=lv; pr.pa=modpar26(3); pr.pa.oae=1; pr.pa.rough_amp=ra;
  S=tdm26('wnr1',pr,0,0);
  od=S.od; tpk=S.tpk; reclen=od.t(end);
  [~,ip]=max(od.env); tpeak=od.t(ip);
  fprintf('%9.1e  %7.1f   %9.2f    %9.2f    %8.2f   %5.1f\n', ...
     ra, S.oam, S.oae, tpeak, od.gd, od.snr);
  if (ra==ras(end)), odbig=od; end
end
fprintf('\nreference: WNR forward latency = %.2f ms  ->  2xWNR = %.2f ms\n', tpk, 2*tpk);
fprintf('record length = %.2f ms (record-center = %.2f ms, the noise-centroid value)\n', reclen, reclen/2);
save('rough_sweep_tboae.mat','odbig','ras','fr','lv','tpk');
disp('ROUGH_SWEEP_DONE');
