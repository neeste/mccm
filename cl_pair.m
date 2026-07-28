% OHC as an INTERNAL FORCE PAIR between BM and RL (anatomically correct).
% Which sign amplifies?  Test both, plus the passive control (ohcgain=0).
% 3-chamber ref @2kHz: lat 6.44/4.72/4.40/2.16, WNR ratios 3.08/2.91/10.86
lv=[20 40 60 80];
runs={ {'passive (gain=0)',0,+1}, {'pair sgn=+1',1,+1}, {'pair sgn=-1',1,-1} };
for c=1:numel(runs)
    nm=runs{c}{1}; g=runs{c}{2}; sg=runs{c}{3};
    fprintf('\n%s\n',nm);
    fprintf('  lvl  lat(ms)  d2/hbmx      max|WNR|  ratio\n');
    prev=NaN;
    for j=1:numel(lv)
        p.fr=2; p.lv=lv(j); p.pa=modpar26(4); p.pa.hbmode='bm';
        p.pa.ohcgain=g; p.pa.ohcsgn=sg;
        try
            S=tdm26('wnr1',p,0,0); d=S.dgn; mx=max(S.wnr);
            r=NaN; if (isfinite(prev)&&prev>0), r=mx/prev; end
            fprintf('  %3.0f %7.2f %8.2f  %11.3e %6.2f\n',lv(j),S.tpk,d.ratio,mx,r);
            prev=mx;
        catch e
            fprintf('  %3.0f  FAILED: %s\n',lv(j),e.message);
        end
    end
end
fprintf('\nAMPLIFYING sign => larger WNR than passive, ratio<10, longer latency.\n');
fprintf('DAMPING sign    => smaller WNR than passive.\n');
disp('CL_PAIR_DONE');
