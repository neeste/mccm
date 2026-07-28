function a2_chsz2_test
% (a) Reduce the 2nd-chamber fluid coupling chsz(2): does the basal ~5.5 ms
% reflection collapse while the CF map (WNR) and the CF peak (~2xWNR) survive?
% Full roughness, 0.5 kHz / 60 dB.
fr=0.5;
fprintf('\n=== 2nd-chamber coupling chsz(2) sweep @ %.1f kHz / 60 dB ===\n',fr);
fprintf('  chsz2   WNR    basal[3-8](t:h)   CF[11-17](t:h)   dominant\n');
for c2=[1 0.5 0.25 0.1]
  pr.fr=fr; pr.lv=60; pr.pa=modpar26(3); pr.pa.oae=1; pr.pa.rough_amp=3e-2;
  pr.pa.chsz=[1 c2 1];
  S=tdm26('wnr1',pr,0,0); od=S.od; t=od.t; env=od.env; g=max(env(t>=1&t<=35));
  [bt,bh]=pkin(t,env,3,8,g); [ct,ch]=pkin(t,env,11,17,g);
  dom='basal'; if(ch>bh),dom='CF'; end
  fprintf('  %5.2f   %5.2f    %5.2f:%.2f     %6.2f:%.2f      %s\n',c2,S.tpk,bt,bh,ct,ch,dom);
end
disp('A2_CHSZ2_DONE');
end
function [tp,hp]=pkin(t,env,lo,hi,g)
w=t>=lo&t<=hi; tt=t(w); ee=env(w); [hp,i]=max(ee); tp=tt(i); hp=hp/g;
end
