% Does the CL chamber affect the TBOAE emission, even though it cannot amplify?
% Matched: rough_amp=3e-2, hbmode='bm', 60 dB.  Emission = ped_rough - ped_smooth.
%   3-chamber            : reference
%   4-chamber ohcgain=0  : CL present, pump OFF  -> passive fluid path only
%   4-chamber ohcgain=1  : CL pump ACTIVE (stable) -> adds an active source
L=load('parfit26_recip.mat');
cfg={ {'3-chamber',       L.R.pa,        []}, ...
      {'4-ch CL pump OFF', modpar26(4),  0 }, ...
      {'4-ch CL pump ON',  modpar26(4),  1 } };
fprintf('\n%-18s %5s | %6s %6s | %8s | %s\n','config','f','WNR','2xWNR','emag(dB)','OAE lat / peak nearest 2xWNR');
fprintf('%s\n',repmat('-',1,96));
for c=1:numel(cfg)
    nm=cfg{c}{1}; pa=cfg{c}{2}; g=cfg{c}{3};
    for fr=[0.5 2]
        p.fr=fr; p.lv=60; p.pa=pa;
        p.pa.hbmode='bm'; p.pa.oae=1; p.pa.rough_amp=3e-2;
        if (~isempty(g)), p.pa.ohcgain=g; p.pa.ohcsgn=+1; end
        try
            S=tdm26('wnr1',p,0,0);
            if (isempty(S.od) || ~isfinite(S.tpk))
                fprintf('%-18s %5.2f | (no valid response)\n',nm,fr); continue; end
            od=S.od; t=od.t; env=od.env; w=t>=1&t<=35; g2=max(env(w));
            ii=find(w); ii=ii(2:end-1);
            lm=ii(env(ii)>env(ii-1)&env(ii)>=env(ii+1));
            tgt=2*S.tpk; s='';
            if (~isempty(lm))
                [~,jn]=min(abs(t(lm)-tgt)); k=lm(jn);
                s=sprintf('%.2f ms | near2x %.2f (%.2f) d=%+.2f', S.oae, t(k), env(k)/g2, t(k)-tgt);
            end
            fprintf('%-18s %5.2f | %6.2f %6.2f | %8.1f | %s\n', nm, fr, S.tpk, tgt, S.oam, s);
        catch e
            fprintf('%-18s %5.2f | FAILED: %s\n', nm, fr, e.message);
        end
    end
end
fprintf('\nCL matters if emag or the emission structure differs from the 3-chamber,\n');
fprintf('and pump ON vs OFF separates an ACTIVE source from a passive fluid path.\n');
disp('CL_OAE_DONE');
