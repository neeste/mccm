function z0slope_sweep
% Can ANY z0slope move dominance from the fixed basal peak (~5.5 ms) to the CF
% peak (~2xWNR ~ 14 ms)?  If yes, the basal peak IS a z0-gradient reflection and
% we just need the right slope; if no slope works, it is something else.
pa=modpar26(3);
d=pa.bwe-pa.ace; sz_def=(pa.k1e-d)/2 - pa.ace;
fprintf('k1e=%.4f ace=%.4f bwe=%.4f  ->  default z0unif sz = %.4f\n',...
        pa.k1e, pa.ace, pa.bwe, sz_def);
fr=0.5;
szs=unique([0, sz_def*[1 2 4 8], -2*sz_def]);
fprintf('\n@ %.1f kHz / 60 dB   basal pk in [3,8]ms   CF pk in [11,17]ms (~2xWNR)\n',fr);
fprintf('  z0slope    WNR    basal(t:h)     CF(t:h)     CF>basal?\n');
for sz=szs
  pr.fr=fr; pr.lv=60; pr.pa=modpar26(3); pr.pa.oae=1; pr.pa.rough_amp=3e-2;
  if (sz~=0), pr.pa.z0unif=1; pr.pa.z0slope=sz; end
  S=tdm26('wnr1',pr,0,0); od=S.od; t=od.t; env=od.env; g=max(env(t>=1&t<=35));
  [bt,bh]=pkin(t,env,3,8,g); [ct,ch]=pkin(t,env,11,17,g);
  fprintf('  %7.3f  %5.2f   %5.2f:%.2f   %6.2f:%.2f      %d\n', sz,S.tpk,bt,bh,ct,ch,ch>bh);
end
disp('Z0SLOPE_SWEEP_DONE');
end
function [tp,hp]=pkin(t,env,lo,hi,g)
w=t>=lo&t<=hi; tt=t(w); ee=env(w); [hp,i]=max(ee); tp=tt(i); hp=hp/g;
end
