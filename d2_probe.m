% Is d2 (bundle displacement) reaching hbmx?  gam compresses only when |d2|>hbmx.
fr=2; lv=[40 80];
L=load('parfit26_recip.mat');
cfg={ {'3-chamber',L.R.pa}, {'4-chamber',modpar26(4)} };
fprintf('\n%-11s %4s   %11s %11s %11s   %s\n','model','lvl','max|d1|','max|d2|','hbmx','d2/hbmx');
for c=1:2
  nm=cfg{c}{1}; pa=cfg{c}{2}; pa.hbmode='bm';
  for j=1:numel(lv)
    p.fr=fr; p.lv=lv(j); p.pa=pa;
    S=tdm26('wnr1',p,0,0); d=S.dgn;
    fprintf('%-11s %4.0f   %11.3e %11.3e %11.3e   %8.2f\n', nm, lv(j), d.d1mx, d.d2mx, d.hbmx, d.ratio);
  end
end
fprintf('\nd2/hbmx > 1 => compression engages;  << 1 => model runs LINEAR\n');
disp('D2_PROBE_DONE');
