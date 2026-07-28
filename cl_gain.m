% With sgn=-1 established as the AMPLIFYING sense (but unstable at full gain),
% find the stable amplifying regime, using the ENERGY diagnostic as the arbiter:
%   ohcP > 0  => OHC force opposes damping and injects energy (true amplifier)
%   ohcP < 0  => dissipating (wrong sign or wrong k5 phase)
% Stable amplification = ohcP>0, WNR above passive, ratio<10, latency responding.
lv=[20 60];
fprintf('\n  sgn  gain |  lvl  lat(ms)  d2/hbmx     max|WNR|      ohcP        ohcW\n');
fprintf('%s\n',repmat('-',1,78));
cfgs={{+1,0},{-1,0.02},{-1,0.05},{-1,0.1},{-1,0.3},{+1,1.0}};
for c=1:numel(cfgs)
    sg=cfgs{c}{1}; g=cfgs{c}{2};
    for j=1:numel(lv)
        p.fr=2; p.lv=lv(j); p.pa=modpar26(4); p.pa.hbmode='bm';
        p.pa.ohcsgn=sg; p.pa.ohcgain=g;
        try
            S=tdm26('wnr1',p,0,0); d=S.dgn;
            fprintf('  %+d %5.2f | %3.0f %7.2f %8.2f  %11.3e %+11.3e %+11.3e\n', ...
                    sg,g,lv(j),S.tpk,d.ratio,max(S.wnr),d.ohcP,d.ohcW);
        catch e
            fprintf('  %+d %5.2f | %3.0f  FAILED: %s\n',sg,g,lv(j),e.message);
        end
    end
end
fprintf('\n(gain=0 rows are the passive control: ohcP must be ~0)\n');
disp('CL_GAIN_DONE');
