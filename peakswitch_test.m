% Is the 40->60 dB latency collapse a DETECTOR PEAK SWITCH (late BM peak losing
% to the early 2nd-DOF shoulder)?  Show the shoulder grid, then the actual WNR
% peak structure at 0.5 kHz for 40 vs 60 dB.
L=load('tmp_surface_check.mat'); m=L.R.mf;
fprintf('shoulder ratio grid (rows f=0.5/1/2/4 kHz, cols 20/40/60/80 dB):\n');
disp(round(m.sho,3));
fprintf('latency grid:\n'); disp(round(m.lat,2));
fprintf('\n--- WNR peak structure, 0.5 kHz ---\n');
for lv=[40 60]
  pr.fr=0.5; pr.lv=lv; pr.pa=L.R.pa;
  S=tdm26('wnr1',pr,0,0);
  w=S.wnr(:); dt=S.dtms; t=(0:numel(w)-1)'*dt;
  win=t>=1 & t<=20; ii=find(win); ii=ii(2:end-1);
  lm=ii(w(ii)>w(ii-1) & w(ii)>=w(ii+1));
  [~,o]=sort(w(lm),'descend'); lm=lm(o);
  fprintf('%3.0f dB: detector latency=%.2f ms | top peaks (ms:rel):', lv, S.tpk);
  for q=1:min(4,numel(lm)), fprintf('  %.2f(%.2f)', t(lm(q)), w(lm(q))/w(lm(1))); end
  fprintf('\n');
end
disp('PEAKSWITCH_DONE');
