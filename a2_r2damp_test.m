function a2_r2damp_test
% (a) Suppress the basal 2nd-DOF over-reflection by damping the 2nd radial DOF
% (r2o). Goal: kill the ~5.5 ms basal peak while (i) preserving the CF map (WNR)
% and (ii) letting the CF coherent reflection (~2xWNR~13.6 ms) emerge dominant.
% Full (physical) roughness, 0.5 kHz / 60 dB.
fr=0.5; p0=modpar26(3); r2o0=p0.r2o;
fprintf('\n=== 2nd-DOF damping (r2o x mult) @ %.1f kHz / 60 dB ===\n',fr);
fprintf('  r2mult   WNR    basal[3-8](t:h)   CF[11-17](t:h)   dominant\n');
for mult=[1 2 4 8]
  pr.fr=fr; pr.lv=60; pr.pa=modpar26(3); pr.pa.oae=1; pr.pa.rough_amp=3e-2;
  pr.pa.r2o=r2o0*mult;
  S=tdm26('wnr1',pr,0,0); od=S.od; t=od.t; env=od.env; g=max(env(t>=1&t<=35));
  [bt,bh]=pkin(t,env,3,8,g); [ct,ch]=pkin(t,env,11,17,g);
  dom='basal'; if(ch>bh),dom='CF'; end
  fprintf('  %5.0f   %5.2f    %5.2f:%.2f     %6.2f:%.2f      %s\n',mult,S.tpk,bt,bh,ct,ch,dom);
end
disp('A2_R2DAMP_DONE');
end
function [tp,hp]=pkin(t,env,lo,hi,g)
w=t>=lo&t<=hi; tt=t(w); ee=env(w); [hp,i]=max(ee); tp=tt(i); hp=hp/g;
end
