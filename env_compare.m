function env_compare
% TEST B: does the 1-CHAMBER also produce a TBOAE coherent-reflection peak at
% ~2x the forward latency, as the 3-chamber does?  Matched conditions: same
% roughness (3e-2), same level (60 dB), BM-peak forward reference (hbmode='bm').
cfgs={'refit_c1broad.mat','parfit26_recip.mat'};
for c=1:numel(cfgs)
    if (~exist(cfgs{c},'file')), fprintf('%s missing\n',cfgs{c}); continue; end
    L=load(cfgs{c}); pa0=L.R.pa; nch=1; if (isfield(pa0,'m')), nch=pa0.m; end
    fprintf('\n=== %s   (nch=%d) ===\n', cfgs{c}, nch);
    fprintf(' f(kHz)   WNR  2xWNR |  dominant peak  |  peak nearest 2xWNR   | emag(dB)\n');
    for fr=[0.5 1 2 4]
        pr.fr=fr; pr.lv=60; pr.pa=pa0;
        pr.pa.oae=1; pr.pa.rough_amp=3e-2; pr.pa.hbmode='bm';
        S=tdm26('wnr1',pr,0,0);
        if (isempty(S.od) || ~isfinite(S.tpk))
            fprintf('%6.2f   (no valid response)\n',fr); continue; end
        od=S.od; t=od.t; env=od.env;
        win=t>=1 & t<=35; g=max(env(win));
        ii=find(win); ii=ii(2:end-1);
        lm=ii(env(ii)>env(ii-1) & env(ii)>=env(ii+1));
        if (isempty(lm)), fprintf('%6.2f   (no envelope peaks)\n',fr); continue; end
        [~,o]=sort(env(lm),'descend'); lmS=lm(o);
        tgt=2*S.tpk;
        [~,jn]=min(abs(t(lm)-tgt)); inear=lm(jn);
        fprintf('%6.2f %5.2f %6.2f | %6.2f (%.2f)  | %6.2f (%.2f) d=%+5.2f | %7.1f\n', ...
            fr, S.tpk, tgt, t(lmS(1)), env(lmS(1))/g, t(inear), env(inear)/g, t(inear)-tgt, S.oam);
    end
end
fprintf('\nA CF coherent reflection shows as a peak at d~0 (i.e. sitting at 2xWNR) with meaningful height.\n');
disp('ENV_COMPARE_DONE');
