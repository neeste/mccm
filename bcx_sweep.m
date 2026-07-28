function bcx_sweep
% Test the basal-boundary-mismatch hypothesis: scale ONLY the stapes boundary
% admittance (2*alfx -> 2*alfx*bcx). 2*alfx/L1_c = 0.018, so matched ~ bcx=57.
% Prediction: as bcx rises, the near-rigid base becomes absorbing -> the basal
% reflection DROPS and the emission magnitude RISES (more couples out to ped).
fr=0.5;
fprintf('\n=== basal-boundary admittance scale (bcx) @ %.1f kHz / 60 dB ===\n',fr);
fprintf('  bcx    WNR    basal[3-8](t:h)   CF[11-17](t:h)  dom    emag(dB)\n');
for bcx=[1 8 30 57 114]
  pr.fr=fr; pr.lv=60; pr.pa=modpar26(3); pr.pa.oae=1; pr.pa.rough_amp=3e-2; pr.pa.bcx=bcx;
  S=tdm26('wnr1',pr,0,0); od=S.od; t=od.t; env=od.env; g=max(env(t>=1&t<=35));
  [bt,bh]=pkin(t,env,3,8,g); [ct,ch]=pkin(t,env,11,17,g);
  dom='bas'; if(ch>bh),dom='CF '; end
  fprintf('  %4.0f   %5.2f    %5.2f:%.2f     %6.2f:%.2f   %s   %6.1f\n',bcx,S.tpk,bt,bh,ct,ch,dom,S.oam);
end
disp('BCX_SWEEP_DONE');
end
function [tp,hp]=pkin(t,env,lo,hi,g)
w=t>=lo&t<=hi; tt=t(w); ee=env(w); [hp,i]=max(ee); tp=tt(i); hp=hp/g;
end
