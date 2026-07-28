function recip_test
% Make the stapes port bi-directional (reverse coupling = transpose of forward,
% both carrying cp.alfx) and test the three things it should fix:
%   (1) emission magnitude emag should RISE (energy now couples outward),
%   (2) the spurious basal reflection should DROP (base properly terminated),
%   (3) the forward response WNR should SURVIVE (drive path unchanged).
fprintf('\n=== bi-directional (reciprocal) stapes coupling ===\n');
fprintf(' f(kHz) recip  WNR    OAE    2xWNR  basal[3-8]  CF-region   emag(dB)\n');
for fr=[0.5 2]
  for rc=[0 1]
    pr.fr=fr; pr.lv=60; pr.pa=modpar26(3); pr.pa.oae=1; pr.pa.rough_amp=3e-2;
    pr.pa.me_recip=rc;   % explicit: 0 = legacy asymmetric port, 1 = reciprocal (now default)
    S=tdm26('wnr1',pr,0,0);
    if (~isfinite(S.tpk) || isempty(S.od))
      fprintf('%6.2f  %d     (no valid response)\n',fr,rc); continue; end
    od=S.od; t=od.t; env=od.env; g=max(env(t>=1&t<=35));
    lo=max(1,2*S.tpk-3); hi=2*S.tpk+3;
    [bt,bh]=pkin(t,env,3,8,g); [ct,ch]=pkin(t,env,lo,hi,g);
    fprintf('%6.2f  %d   %5.2f  %5.2f  %5.2f  %5.2f:%.2f  %5.2f:%.2f  %7.1f\n', ...
            fr, rc, S.tpk, S.oae, 2*S.tpk, bt,bh, ct,ch, S.oam);
  end
end
disp('RECIP_TEST_DONE');
end
function [tp,hp]=pkin(t,env,lo,hi,g)
w=t>=lo&t<=hi; if(~any(w)), tp=NaN; hp=NaN; return; end
tt=t(w); ee=env(w); [hp,i]=max(ee); tp=tt(i); hp=hp/g;
end
