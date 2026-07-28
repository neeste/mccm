% PHASE 2: does routing OHC somatic motility into the cortilymph pump (d3)
% smooth the COMPRESSION CLIFF?  The cliff is the per-20 dB latency step: every
% 3-chamber configuration and every transduction drive gave ~53-55%, against the
% ~28%/20 dB the data requires.  Test the level series at one frequency for
% several split fractions.
%   ohcsplit=0 must reproduce validated Phase 1 (regression gate: 1.90 ms @2k/60)
fr=2; lv=[20 40 60 80]; fsps=[0 0.5 1.0];
fprintf('\nPHASE 2 -- OHC somatic split into the cortilymph pump (%g kHz)\n',fr);
fprintf('%8s |%s | %s\n','ohcsplit',sprintf('%8.0fdB',lv),'per-20dB steps');
fprintf('%s\n',repmat('-',1,72));
for fsp=fsps
    lat=nan(1,numel(lv));
    for j=1:numel(lv)
        p.fr=fr; p.lv=lv(j); p.pa=modpar26(4); p.pa.hbmode='bm'; p.pa.ohcsplit=fsp;
        try
            S=tdm26('wnr1',p,0,0); lat(j)=S.tpk;
        catch e
            fprintf('  fsp=%.2f lv=%d FAILED: %s\n',fsp,lv(j),e.message);
        end
    end
    s='';
    for j=1:numel(lv)-1
        if (isfinite(lat(j))&&isfinite(lat(j+1)))
            s=[s sprintf('%7.1f%%',100*(lat(j+1)-lat(j))/lat(j))];
        else
            s=[s '      NaN'];
        end
    end
    fprintf('%8.2f |%s | %s\n', fsp, sprintf('%10.2f',lat), s);
end
fprintf('\ntarget: uniform ~-28%%/20dB.  3-chamber and all drives gave ~-53 to -55%% worst step.\n');
fprintf('ohcsplit=0 row must match Phase 1 (1.90 ms at 60 dB).\n');
disp('CL_PHASE2_DONE');
