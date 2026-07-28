function rough_window_test
% If the fixed ~5 ms basal peak is BASAL-ROUGHNESS reflection, restricting the
% roughness apically (rough_xlo = basal cutoff fraction) should make it VANISH
% while the CF peak (~2xWNR ~ 14 ms at 0.5 kHz) survives.
fr=0.5;
xlos=[0 0.2 0.4 0.6];   % 0 = full-length roughness (baseline)
fprintf('\n=== roughness basal-cutoff @ %.1f kHz / 60 dB ===\n',fr);
fprintf('  rough_xlo   WNR    basal[3-8](t:h)   CF[11-17](t:h)   dominant\n');
for xlo=xlos
  pr.fr=fr; pr.lv=60; pr.pa=modpar26(3); pr.pa.oae=1; pr.pa.rough_amp=3e-2;
  if (xlo>0), pr.pa.rough_xlo=xlo; end
  S=tdm26('wnr1',pr,0,0); od=S.od; t=od.t; env=od.env; g=max(env(t>=1&t<=35));
  [bt,bh]=pkin(t,env,3,8,g); [ct,ch]=pkin(t,env,11,17,g);
  dom='basal'; if (ch>bh), dom='CF'; end
  fprintf('  %6.2f    %5.2f    %5.2f:%.2f     %6.2f:%.2f      %s\n', xlo,S.tpk,bt,bh,ct,ch,dom);
end
disp('ROUGH_WINDOW_DONE');
end
function [tp,hp]=pkin(t,env,lo,hi,g)
w=t>=lo&t<=hi; tt=t(w); ee=env(w); [hp,i]=max(ee); tp=tt(i); hp=hp/g;
end
