% Energy diagnostic (robust). Does the OHC force inject energy, and is there a
% STABLE amplifying regime?  ohcP>0 = injects. gain=0 is the validity control.
lv=[20 60];
fprintf('\n  sgn  gain |  lvl  lat(ms)  d2/hbmx     max|WNR|        ohcP   nonfin\n');
fprintf('%s\n',repmat('-',1,76));
cfgs={{+1,0},{-1,0.1},{-1,0.3},{-1,0.6},{-1,1.0},{+1,1.0}};
for c=1:numel(cfgs)
    sg=cfgs{c}{1}; g=cfgs{c}{2};
    for j=1:numel(lv)
        p.fr=2; p.lv=lv(j); p.pa=modpar26(4); p.pa.hbmode='bm';
        p.pa.ohcsgn=sg; p.pa.ohcgain=g;
        try
            S=tdm26('wnr1',p,0,0); d=S.dgn;
            fprintf('  %+d %5.2f | %3.0f %7.2f %8.2f  %11.3e %+11.3e %6d\n', ...
                    sg,g,lv(j),S.tpk,d.ratio,max(S.wnr),d.ohcP,d.ohcNaN);
        catch e
            fprintf('  %+d %5.2f | %3.0f FAILED: %s\n',sg,g,lv(j),e.message);
        end
    end
end
fprintf('\ngain=0 control: ohcP must be ~0.  ohcP>0 => energy injected.\n');
disp('CL_GAIN2_DONE');
