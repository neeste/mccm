function bmpeak_test
% Two questions:
% (1) CORRECTNESS: is hbmode='bm' a PURE measurement switch?  The eardrum level
%     mlv is a mechanical observable -- it must be bit-identical between the full
%     drive (v1-a*v2) and the BM-only drive (v1) if the mechanics are untouched.
% (2) DOES IT FIX THE DETECTOR: is the BM-only WNR single-peaked, with a SMOOTH
%     level dependence (no 40->60 dB collapse)?
L=load('parfit26_recip.mat'); pa0=L.R.pa;
for fr=[0.5 2]
  fprintf('\n=== %.1f kHz ===\n',fr);
  fprintf(' lvl   mlv(full)  mlv(bm)    |dMLV|    lat(full)  lat(bm)   BM-WNR peaks (ms:rel)\n');
  for lv=[20 40 60 80]
    pr.fr=fr; pr.lv=lv; pr.pa=pa0;
    Sf=tdm26('wnr1',pr,0,0);
    pr2.fr=fr; pr2.lv=lv; pr2.pa=pa0; pr2.pa.hbmode='bm';
    Sb=tdm26('wnr1',pr2,0,0);
    w=Sb.wnr(:); dt=Sb.dtms; t=(0:numel(w)-1)'*dt;
    win=t>=1&t<=20; ii=find(win); ii=ii(2:end-1);
    lm=ii(w(ii)>w(ii-1)&w(ii)>=w(ii+1)); [~,o]=sort(w(lm),'descend'); lm=lm(o);
    s='';
    for q=1:min(3,numel(lm)), s=[s sprintf('  %.2f(%.2f)',t(lm(q)),w(lm(q))/w(lm(1)))]; end
    fprintf(' %3.0f   %8.2f  %8.2f  %9.2e  %8.2f  %7.2f  %s\n', ...
            lv, Sf.mlv, Sb.mlv, abs(Sf.mlv-Sb.mlv), Sf.tpk, Sb.tpk, s);
  end
end
disp('BMPEAK_DONE');
