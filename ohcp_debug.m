p.fr=2; p.lv=60; p.pa=modpar26(4); p.pa.hbmode='bm'; p.pa.ohcsgn=-1; p.pa.ohcgain=0.3;
S=tdm26('wnr1',p,0,0); d=S.dgn;
fprintf('non-finite ohcp samples: %d of %d\n', d.ohcNaN, d.ohcN);
fprintf('ohcP=%+.4e   ohcW=%+.4e\n', d.ohcP, d.ohcW);
fprintf('lat=%.2f  d2/hbmx=%.2f  maxWNR=%.3e\n', S.tpk, d.ratio, max(S.wnr));
disp('OHCP_DEBUG_DONE');
