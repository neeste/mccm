% Inspect the emission-envelope peak structure: is there a fixed basal peak
% (~5 ms, freq-independent) AND a masked CF peak (~2x forward, freq-dependent)?
for fr=[0.5 1 2 4]
  pr.fr=fr; pr.lv=60; pr.pa=modpar26(3); pr.pa.oae=1; pr.pa.rough_amp=3e-2;
  S=tdm26('wnr1',pr,0,0); od=S.od; env=od.env; t=od.t;
  w=t>=1 & t<=35; ii=find(w); ii=ii(2:end-1);
  lm=ii(env(ii)>env(ii-1) & env(ii)>=env(ii+1));
  [~,ord]=sort(env(lm),'descend'); lm=lm(ord);
  emx=max(env(w));
  fprintf('\n%.1f kHz/60: WNR=%.2f  2xWNR=%.2f  chosenOAE=%.2f   top env peaks (ms : rel-height):\n',...
          fr, S.tpk, 2*S.tpk, S.oae);
  for q=1:min(5,numel(lm)), fprintf('     %5.2f ms : %.2f\n', t(lm(q)), env(lm(q))/emx); end
end
disp('ENV_DUMP_DONE');
