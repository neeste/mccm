% Is the 4-chamber running LINEAR (compression never engaging)?
% If WNR amplitude scales by exactly x10 per 20 dB, the response is linear and
% the OHC nonlinearity (gam via hbt=max(|d2|/hbmx,1)) is not being triggered.
% The 3-chamber is the control: it must show clear compression (<x10).
fr=2; lv=[20 40 60 80];
L=load('parfit26_recip.mat');
cfg={ {'3-chamber', L.R.pa}, {'4-chamber', modpar26(4)} };
for c=1:numel(cfg)
    nm=cfg{c}{1}; pa=cfg{c}{2}; pa.hbmode='bm';
    fprintf('\n%s @ %g kHz\n', nm, fr);
    fprintf('  lvl    lat(ms)      max|WNR|     ratio/20dB   (linear = 10.00)\n');
    prev=NaN;
    for j=1:numel(lv)
        p.fr=fr; p.lv=lv(j); p.pa=pa;
        try
            S=tdm26('wnr1',p,0,0); mx=max(S.wnr);
            r=NaN; if (isfinite(prev)&&prev>0), r=mx/prev; end
            fprintf('  %3.0f   %7.2f   %11.3e   %10.2f\n', lv(j), S.tpk, mx, r);
            prev=mx;
        catch e
            fprintf('  %3.0f   FAILED: %s\n', lv(j), e.message);
        end
    end
end
fprintf('\nratio ~10.00 at every step => LINEAR (no compression engaging)\n');
fprintf('ratio < 10 => compression active (3-chamber should show this)\n');
disp('CL_LINEAR_DONE');
