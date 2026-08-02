% C3PROBE -- do ace and hbsc carry to the 3-chamber model?
%
% Everything established 2026-08-01 -- ace as the slope lever, hbsc as the
% compression slope that puts level_c in band -- was measured at nch=1 ONLY. The
% multi-chamber model is the capstone's MECHANISM model, so before spending the
% Banff window fitting it, check the two levers actually do the same thing there.
%
% THE TRAP THIS AVOIDS. parfit26 drops a warm start whose chamber count does not
% match: `if (numel(pa0.chsz) ~= numel(base.chsz)), pa0=base; end`
% (parfit26.m:157). The nch=1 result has 2 chsz entries, modpar26(3) has 3, so
% handing the nch=1 pa to parfit26(3) would SILENTLY fall back to the baseline
% and throw away ace and hbsc without a word. That already bit this project once
% -- sweep1b passed no 'warm', hit the same check, and ran from the baseline
% while appearing to warm-start. So: start from the best nch=3 fit on record and
% transplant ONLY the two levers onto it.
%
% parfit26_nch3.mat is the best nch=3 warm start (maperr 147.1, maxRe_osc -220.3,
% and sub-critical -- all six nch=3 starts on disk were checked 2026-08-01).
%
% Row 1 is the ANCHOR: untouched, must reproduce maperr ~147.1 and osc ~-220.7.
% If it does not, the file is not what the record says and nothing below counts.
t0=tic;
L1=load('/Users/neely/mccm_runs/longfit_seg2.mat'); p1=L1.R.pa;
L3=load('parfit26_nch3.mat');  p3=L3.R.pa;
ACE=p1.ace; HBSC=p1.hbsc;
fprintf('\n== do ace/hbsc carry to nch=3? ==\n');
fprintf('   from nch=1 seg2: ace %+.4f  hbsc %.4f  (there: slope 0.4130, lvl_c 3.816, maperr 86.41)\n', ACE, HBSC);
fprintf('   nch=3 base parfit26_nch3.mat: ace %+.4f  hbsc %.4f  chsz %s\n\n', ...
        p3.ace, p3.hbsc, mat2str(round(p3.chsz,4)));
fprintf('   config                    slope   lvl_c   maperr    amp     osc    |  Jslp   Jlcb   Jmap  Jgain  Josc =    J\n');
fprintf('   ----------------------------------------------------------------------------------------------------------\n');
cfg = { p3.ace, p3.hbsc, 'ANCHOR nch=3 as-is';
        p3.ace, HBSC,    '+ hbsc only';
        ACE,    p3.hbsc, '+ ace only';
        ACE,    HBSC,    '+ both' };
MAPTOL = 150;   % the nch=3 default; nch=1 used 105
best=Inf; bi=NaN; R=struct();
for j=1:size(cfg,1)
    pa=p3; pa.ace=cfg{j,1}; pa.hbsc=cfg{j,2};
    sl=NaN; lc=NaN; mp=NaN; ag=NaN; mo=NaN;
    try, m=abr_metric(pa,false); if (m.ok), sl=m.slope; lc=m.level_c; end, catch, end
    try
        S=score26(pa,'fast',false); mp=S.maperr; mo=S.maxRe_osc;
        if (isfield(S,'amp_gain')&&~isempty(S.amp_gain)), ag=double(S.amp_gain(1)); end
    catch, end
    Js=abs(sl-0.413); Jl=0.10*(max(0,3.5-lc)+max(0,lc-6.5));
    Jp=0.001*max(0,mp-MAPTOL); Jg=0.01*max(0,40-ag); Jo=0.005*max(0,mo+40);
    J=Js+Jl+Jp+Jg+Jo;
    if (isfinite(J) && J<best), best=J; bi=j; end
    fprintf('   %-22s  %7.4f  %6.3f  %7.2f  %+6.2f  %7.1f | %6.4f %6.4f %6.4f %6.4f %6.4f = %6.4f\n', ...
            cfg{j,3}, sl, lc, mp, ag, mo, Js, Jl, Jp, Jg, Jo, J);
    R(j).ace=cfg{j,1}; R(j).hbsc=cfg{j,2}; R(j).slope=sl; R(j).level_c=lc;
    R(j).maperr=mp; R(j).amp=ag; R(j).osc=mo; R(j).J=J;
end
save('c3probe.mat','R','cfg','ACE','HBSC');
fprintf('\n   best row %d, J %.4f   |   %.1f min\n', bi, best, toc(t0)/60);
fprintf('   READ: does the slope fall toward 0.413 and level_c rise toward 3.5,\n');
fprintf('         as they did at nch=1? And does osc stay well below -40?\n');
disp('C3PROBE_DONE');
