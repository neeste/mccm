% HBSWEEP2 -- where does hbsc saturate, and does hbmx UP help?
%
% hbsweep found the first non-trading lever of the session: raising hbsc from
% 0.04 to 0.32 at ace=+0.40 lifted level_c 2.818 -> 3.253 AND lowered the slope
% 0.5067 -> 0.4827, with maperr pinned at 90.52 on every row. Jslp+Jlcb went
% 0.1619 -> 0.0944, against 0.1928 at the fitted optimum.
%
% Two open ends, both cheap:
%   1. level_c gained only 0.010 between hbsc 0.16 and 0.32 while the slope kept
%      falling. That looks like level_c SATURATING while the slope keeps going --
%      if so the two terms stop moving together and the useful range ends. Push
%      hbsc to 0.64 and 1.28 to find the knee rather than assume it.
%   2. hbmx DOWN was firmly wrong (level_c 1.700 at 6e-10, 1.300 at 6e-11), which
%      makes hbmx UP the untested direction. d2mx/hbmx ~ 21 at baseline, so
%      raising hbmx moves the compression knee TOWARD the operating point.
%
% maperr stays the control: fdm26 never sees hbsc or hbmx, so any row where it
% leaves 90.52 invalidates that row.
t0=tic;
L=load('sweep_nch1b.mat'); pa0=L.R.pa; pa0.ace=0.40;
fprintf('\n== hbsc saturation and hbmx UP, at ace=+0.40 ==\n');
fprintf('   best so far: hbsc 0.320 hbmx 6e-09 -> level_c 3.253 slope 0.4827 Jsum 0.0944\n');
fprintf('   band [3.5 6.5]; level_c must reach 3.5 for Jlcb to hit zero.\n\n');
fprintf('    hbsc      hbmx     slope   level_c   Jslp    Jlcb    Jsum  |  maperr(control)\n');
fprintf('   --------------------------------------------------------------------------------\n');
cfg = { 0.320, 6e-9;      % anchor: must reproduce 0.4827 / 3.253
        0.640, 6e-9;
        1.280, 6e-9;
        0.040, 6e-8;
        0.320, 6e-8;
        0.320, 6e-7;
        0.640, 6e-8 };
best=Inf; bj=NaN; R=struct();
for j=1:size(cfg,1)
    pa=pa0; pa.hbsc=cfg{j,1}; pa.hbmx=cfg{j,2};
    sl=NaN; lc=NaN; mp=NaN;
    try, m=abr_metric(pa,false); if (m.ok), sl=m.slope; lc=m.level_c; end, catch, end
    try, S=score26(pa,'fast',false); mp=S.maperr; catch, end
    Js=abs(sl-0.413); Jl=0.10*(max(0,3.5-lc)+max(0,lc-6.5)); Jt=Js+Jl;
    if (isfinite(Jt) && Jt<best), best=Jt; bj=j; end
    fprintf('   %6.3f  %8.2g  %7.4f  %7.3f  %6.4f  %6.4f  %6.4f  |  %7.2f\n', ...
            cfg{j,1}, cfg{j,2}, sl, lc, Js, Jl, Jt, mp);
    R(j).hbsc=cfg{j,1}; R(j).hbmx=cfg{j,2}; R(j).slope=sl; R(j).level_c=lc;
    R(j).maperr=mp; R(j).J=Jt;
end
save('hbsweep2.mat','R','cfg');
fprintf('\n   best Jslp+Jlcb: row %d, %.4f\n', bj, best);
fprintf('   %.1f min.  maperr must be 90.52 on every row.\n', toc(t0)/60);
disp('HBSWEEP2_DONE');
