% AMPCHECK -- is hbsc=0.64 buying the ABR fit by gutting the amplifier?
%
% THE WORRY. hbsc=0.64 at ace=+0.40 gives slope 0.4160 (target 0.413) and
% level_c 3.732 (in band for the first time in this project), J 0.0030 against
% 0.1928. But amp_gain read +44.62 IDENTICALLY at hbsc 0.04 / 0.32 / 0.64 / 1.28,
% and score26.m:56 records a PRIOR bug with exactly that signature -- amp_gain
% "returned +2.46 dB identically for ohcgain 1, 2 and 4 because every row was
% measured at 1". An unmoving guard is the thing to distrust, especially this
% guard, which exists because the model was once caught scoring its best on a
% config with amp_gain -0.01 dB.
%
% WHY IT DID NOT MOVE, found by reading rather than guessing: score26's
% local_tip runs tdm26(0,pa,0,0), a CLICK, and tdm26.m:503 takes pa.hbnl from
% the passed struct. Every sweep row left pa.hbnl at the saved fit's 0, so the
% click was LINEAR while the ABR path ran compressive (tbabr_condition forces
% hbnl=1 at tdm26.m:1536). So amp_gain measured the small-signal amplifier, which
% hbsc cannot affect by construction: below hbmx, log(max(|dh|/hbmx,1)) = 0 and
% gam = cp.gm exactly, whatever hbsc is.
%
% That is the benign reading and it is probably right, but it is still an
% ARGUMENT. Measure instead: force hbnl=1 so the click compresses too, and see
% what the amplifier is worth at a level where compression is active. Rows with
% hbnl=0 are the controls and must reproduce +44.62.
%
% WHAT WOULD BE DAMNING: amp_gain collapsing toward 0 as hbsc rises, i.e. the ABR
% fit bought by switching the amplifier off. WHAT IS FINE: a moderate drop --
% compression REDUCING gain at high level is what compression IS, and the guard's
% 40 dB floor is calibrated on the linear measure.
t0=tic;
L=load('sweep_nch1b.mat'); pa0=L.R.pa; pa0.ace=0.40;
fprintf('\n== amplifier gain with compression actually engaged (ace=+0.40) ==\n');
fprintf('   hbnl=0 rows are CONTROLS: must read +44.62 (the linear amplifier).\n');
fprintf('   gainmin 40 dB. d2mx/hbmx at this operating point was ~21.\n\n');
fprintf('    hbsc   hbnl   amp_gain   amp_d2    contrast   maperr  | reading\n');
fprintf('   ------------------------------------------------------------------\n');
cfg = {0.040,0; 0.640,0; 0.040,1; 0.320,1; 0.640,1; 1.280,1};
R=struct();
for j=1:size(cfg,1)
    pa=pa0; pa.hbsc=cfg{j,1}; pa.hbnl=cfg{j,2};
    ag=NaN; a2=NaN; ct=NaN; mp=NaN;
    try
        S=score26(pa,'fast',false);
        if (isfield(S,'amp_gain')&&~isempty(S.amp_gain)), ag=double(S.amp_gain(1)); end
        if (isfield(S,'amp_d2')&&~isempty(S.amp_d2)), a2=double(S.amp_d2(1)); end
        if (isfield(S,'contrast')), ct=S.contrast; end
        mp=S.maperr;
    catch e
        fprintf('   %5.3f     %d   FAILED: %s\n', cfg{j,1}, cfg{j,2}, e.message(1:min(60,end)));
        continue;
    end
    rd='';
    if (cfg{j,2}==0), rd='control';
    elseif (isfinite(ag))
        if (ag < 1), rd='AMPLIFIER GONE';
        elseif (ag < 40), rd='below gainmin';
        else, rd='amplifier intact'; end
    end
    fprintf('   %5.3f     %d   %+7.2f  %+7.2f   %8.2f  %7.2f  | %s\n', ...
            cfg{j,1}, cfg{j,2}, ag, a2, ct, mp, rd);
    R(j).hbsc=cfg{j,1}; R(j).hbnl=cfg{j,2}; R(j).amp=ag; R(j).amp_d2=a2;
    R(j).contrast=ct; R(j).maperr=mp;
end
save('ampcheck.mat','R','cfg');
fprintf('\n   %.1f min.\n', toc(t0)/60);
disp('AMPCHECK_DONE');
