% After reverting d1 to ST<->SV (CL now a SIDE compartment pumped by d3):
% do the partition displacements recover, does d2 clear hbmx, and does the
% compression engage (WNR ratio < 10 per 20 dB)?
fr=2; lv=[20 40 60 80];
L=load('parfit26_recip.mat');
cfg={{'3-chamber (ref)',L.R.pa},{'4-chamber (CL side)',modpar26(4)}};
for c=1:numel(cfg)
    nm=cfg{c}{1}; pa=cfg{c}{2}; pa.hbmode='bm';
    fprintf('\n%s  @ %g kHz\n',nm,fr);
    fprintf('  lvl  lat(ms)     max|d1|     max|d2|  d2/hbmx      max|WNR|  ratio\n');
    prev=NaN;
    for j=1:numel(lv)
        p.fr=fr; p.lv=lv(j); p.pa=pa;
        try
            S=tdm26('wnr1',p,0,0); d=S.dgn; mx=max(S.wnr);
            r=NaN; if (isfinite(prev)&&prev>0), r=mx/prev; end
            fprintf('  %3.0f %7.2f  %11.3e %11.3e %8.2f  %11.3e %6.2f\n', ...
                    lv(j), S.tpk, d.d1mx, d.d2mx, d.ratio, mx, r);
            prev=mx;
        catch e
            fprintf('  %3.0f  FAILED: %s\n', lv(j), e.message);
        end
    end
end
fprintf('\nd2/hbmx>1 => compression engages;  WNR ratio<10 => nonlinear (10.00 = linear)\n');
disp('CL_REVERT_DONE');
