% Does z0unif (a) preserve the CF map (WNR forward latency unchanged) and
% (b) suppress the fixed ~5 ms basal peak so the CF (2xWNR) peak dominates?
function z0unif_test
for fr=[0.5 2 4]
  a=run1(fr,false);   % no compensation
  b=run1(fr,true);    % z0unif on (default k1-scalar slope)
  fprintf('\n=== %.1f kHz / 60 dB ===  WNR: off=%.2f  on=%.2f ms  (should match: CF map preserved)\n',...
          fr, a.wnr, b.wnr);
  fprintf('   2xWNR = %.2f ms\n', 2*b.wnr);
  showpk('  OFF ', a);
  showpk('  ON  ', b);
end
end

function r=run1(fr,z0)
pr.fr=fr; pr.lv=60; pr.pa=modpar26(3); pr.pa.oae=1; pr.pa.rough_amp=3e-2;
if (z0), pr.pa.z0unif=1; end
S=tdm26('wnr1',pr,0,0);
r.wnr=S.tpk; r.oae=S.oae; r.od=S.od;
end

function showpk(tag,r)
t=r.od.t; env=r.od.env; w=t>=1 & t<=35; ii=find(w); ii=ii(2:end-1);
lm=ii(env(ii)>env(ii-1) & env(ii)>=env(ii+1)); [~,o]=sort(env(lm),'descend'); lm=lm(o);
emx=max(env(w)); s=sprintf('%s chosenOAE=%.2f  peaks:',tag,r.oae);
for q=1:min(4,numel(lm)), s=[s sprintf('  %.2f(%.2f)',t(lm(q)),env(lm(q))/emx)]; end
disp(s);
end
