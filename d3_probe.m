fprintf('\n  sgn  gain |  max|d1|     max|d2|     max|d3|   d3 finite?  maxWNR\n');
for cfg={{+1,0},{-1,0.3},{-1,1.0},{+1,1.0}}
  c=cfg{1}; sg=c{1}; g=c{2};
  p.fr=2; p.lv=60; p.pa=modpar26(4); p.pa.hbmode='bm'; p.pa.ohcsgn=sg; p.pa.ohcgain=g;
  try
    S=tdm26('wnr1',p,0,0); d=S.dgn;
    fprintf('  %+d %5.2f | %10.3e %10.3e %10.3e %8d   %10.3e\n', sg,g,d.d1mx,d.d2mx,d.d3mx,d.d3fin,max(S.wnr));
  catch e, fprintf('  %+d %5.2f | FAILED %s\n',sg,g,e.message); end
end
fprintf('\nd3 finite=0 => the OC-height DOF DIVERGES (invisible to WNR under hbmode=bm)\n');
disp('D3_PROBE_DONE');
