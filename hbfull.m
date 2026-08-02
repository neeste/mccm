% HBFULL -- the hbsc candidates scored on ALL FIVE terms, with amp_gain and osc.
%
% hbsweep2 found hbsc=0.64 at ace=+0.40 giving slope 0.4160 (target 0.413) and
% level_c 3.732 (inside [3.5 6.5] for the first time in this project), for
% Jslp+Jlcb = 0.0030 against 0.1928 at the fitted optimum. Those sweeps measured
% only two of the five terms.
%
% THE TERM I MOST EXPECT TO BREAK IS THE ONE I DID NOT MEASURE. hbsc COMPRESSES
% the amplifier -- gam = cp.gm/(1 + hbsc*log(...)) -- so raising it 16x is
% precisely the thing that should trip `wgain*max(0, 40 - amp_gain)`. Reporting a
% total J without it would repeat this session's mistake twice over.
%
% maxRe_osc should be INERT here: coupeig pins pa.hbnl=0 itself (tdm26.m:1615),
% so compression parameters cannot reach it and osc should equal the ace=+0.40
% value of -83.5 on every row. That is a control, not a result -- if osc moves,
% my model of what coupeig sees is wrong.
%
% THE FORMULA IS THE ONE I GOT WRONG EARLIER, so it is anchored twice here:
%   row 1 is the FITTED OPTIMUM, where parfit26 itself reported J = 0.1928
%   row 2 is ace=+1.00, where the smoke test's own verbterm printed
%          Jm 0.1541 / osc 0 / stat 0 / map 0.0027 / gain 0.0116 = 0.1684
% If both reproduce, Jm = Jslp + Jlcb with no hidden sixth term (wanchor,
% wshoulder, wlevel, wshape all default to 0) and the remaining rows are
% trustworthy. If either misses, nothing below it means anything.
t0=tic;
L=load('sweep_nch1b.mat'); pa0=L.R.pa;
% ace, hbsc, hbmx, label
cfg = { pa0.ace, 0.040, 6e-9,  'ANCHOR fitted optimum (expect J 0.1928)';
        1.000,   0.040, 6e-9,  'ANCHOR smoke point     (expect J 0.1684)';
        0.400,   0.040, 6e-9,  'ace best, stock hbsc';
        0.400,   0.320, 6e-9,  '';
        0.400,   0.640, 6e-9,  'hbsweep2 best';
        0.400,   1.280, 6e-9,  '';
        0.400,   0.320, 6e-8,  '' };
fprintf('\n== full five-term score of the compression candidates ==\n');
fprintf('   maptol 105, gainmin 40, wgain 0.01, wmap 0.001, wcrit 0.005, wlcb 0.10\n\n');
fprintf('    ace   hbsc     hbmx    slope   lvl_c   maperr    amp     osc  |  Jslp   Jlcb   Jmap  Jgain  Josc =    J\n');
fprintf('   ------------------------------------------------------------------------------------------------------\n');
best=Inf; bi=NaN; R=struct();
for j=1:size(cfg,1)
    pa=pa0; pa.ace=cfg{j,1}; pa.hbsc=cfg{j,2}; pa.hbmx=cfg{j,3};
    sl=NaN; lc=NaN; mp=NaN; ag=NaN; mo=NaN;
    try, m=abr_metric(pa,false); if (m.ok), sl=m.slope; lc=m.level_c; end, catch, end
    try
        S=score26(pa,'fast',false); mp=S.maperr; mo=S.maxRe_osc;
        if (isfield(S,'amp_gain')&&~isempty(S.amp_gain)), ag=double(S.amp_gain(1)); end
    catch, end
    Js=abs(sl-0.413);
    Jl=0.10*(max(0,3.5-lc)+max(0,lc-6.5));
    Jp=0.001*max(0,mp-105); Jg=0.01*max(0,40-ag); Jo=0.005*max(0,mo+40);
    J=Js+Jl+Jp+Jg+Jo;
    if (isfinite(J) && j>2 && J<best), best=J; bi=j; end
    fprintf('   %+5.2f  %5.3f  %7.2g  %7.4f  %6.3f  %7.2f  %+6.2f  %7.1f | %6.4f %6.4f %6.4f %6.4f %6.4f = %6.4f  %s\n', ...
            cfg{j,1}, cfg{j,2}, cfg{j,3}, sl, lc, mp, ag, mo, Js, Jl, Jp, Jg, Jo, J, cfg{j,4});
    R(j).ace=cfg{j,1}; R(j).hbsc=cfg{j,2}; R(j).hbmx=cfg{j,3};
    R(j).slope=sl; R(j).level_c=lc; R(j).maperr=mp; R(j).amp=ag; R(j).osc=mo; R(j).J=J;
end
save('hbfull.mat','R','cfg','best','bi');
fprintf('\n   best non-anchor row: %d, J %.4f   (fitted optimum 0.1928)\n', bi, best);
fprintf('   %.1f min.  CHECK THE ANCHORS FIRST: 0.1928 and 0.1684.\n', toc(t0)/60);
fprintf('   CONTROLS: maperr 90.52 (fdm26 blind to hbsc) and osc -83.5 (coupeig pins hbnl=0).\n');
disp('HBFULL_DONE');
