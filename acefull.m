% ACEFULL -- the ace axis scored with the COMPLETE objective, level_c included.
%
% CORRECTING A REAL ERROR. gpsweep, acecheck, acefar and slopeatlas all scored J
% by hand as slope + map + gain + osc. parfit26's Jm also carries a level_c BAND
% term, wlcb=0.10 by default (parfit26.m:153), over lc=[3.5 6.5]:
%
%     Jlcb = wlcb * ( max(0, lclo - level_c) + max(0, level_c - lchi) )
%
% The omission hid itself perfectly at the starting point: the fitted optimum has
% level_c 3.57, INSIDE the band, so Jlcb was exactly 0 and my hand-rolled total
% reproduced parfit26's reported J=0.1928 to four decimals. It only broke once
% ace moved level_c out of the band -- at ace=1.00 the smoke test's own term
% printout reads Jm 0.1541 against my predicted 0.0285, the 0.126 gap being
% 0.10*(3.5-2.24) exactly.
%
% This is the same failure the project keeps hitting: two places holding one
% fact, agreeing only because a quantity happened to be inert. Here the second
% copy was MINE, and the inert quantity was a band term that was satisfied at the
% baseline. The fix is not a better hand formula -- it is to stop keeping a
% second copy. Every row below is scored by calling parfit26's own objective.
%
% WHAT THIS DECIDES: the warm start for an 18-hour run. ace=1.00 was the best
% point on the INCOMPLETE score; with level_c counted the minimum may sit lower,
% since level_c falls monotonically as ace rises and the band penalty grows.
t0=tic;
L=load('sweep_nch1b.mat'); pa0=L.R.pa;
vv=[pa0.ace 0 0.40 0.73 1.00 1.35 1.70];
fprintf('\n== ace axis, FULL objective (level_c band included) ==\n');
fprintf('   lc band [3.5 6.5], wlcb 0.10.  Fitted optimum: J 0.1928, level_c 3.57\n\n');
fprintf('     ace     slope   level_c   maperr    amp     osc   |  Jslp    Jlcb    Jmap   Jgain   Josc  =    J\n');
fprintf('   -------------------------------------------------------------------------------------------------\n');
best=Inf; bace=NaN; R=struct('ace',num2cell(vv));
for j=1:numel(vv)
    pa=pa0; pa.ace=vv(j);
    sl=NaN; lc=NaN; mp=NaN; ag=NaN; mo=NaN;
    try, m=abr_metric(pa,false); if (m.ok), sl=m.slope; lc=m.level_c; end, catch, end
    try
        S=score26(pa,'fast',false); mp=S.maperr; mo=S.maxRe_osc;
        if (isfield(S,'amp_gain')&&~isempty(S.amp_gain)), ag=double(S.amp_gain(1)); end
    catch, end
    Js=abs(sl-0.413);
    Jl=0.10*(max(0,3.5-lc) + max(0,lc-6.5));
    Jp=0.001*max(0,mp-105); Jg=0.01*max(0,40-ag); Jo=0.005*max(0,mo+40);
    J=Js+Jl+Jp+Jg+Jo;
    if (isfinite(J) && J<best), best=J; bace=vv(j); end
    fprintf('   %+7.4f  %7.4f  %7.3f  %7.2f  %+6.2f  %7.1f | %6.4f %6.4f %6.4f %6.4f %6.4f = %6.4f\n', ...
            vv(j), sl, lc, mp, ag, mo, Js, Jl, Jp, Jg, Jo, J);
    R(j).slope=sl; R(j).level_c=lc; R(j).maperr=mp; R(j).amp=ag; R(j).osc=mo; R(j).J=J;
end
save('acefull.mat','R','vv','best','bace');
fprintf('\n   best: ace %+.4f  J %.4f   vs fitted optimum J 0.1928\n', bace, best);
fprintf('   %.1f min.  THIS is the warm start for longfit, not the acefar number.\n', toc(t0)/60);
disp('ACEFULL_DONE');
