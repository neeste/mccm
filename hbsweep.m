% HBSWEEP -- can compression STRENGTH/THRESHOLD restore level_c, the term that
% flattens the ace axis?
%
% WHERE THIS COMES FROM. ace reaches the ABR slope target but drags level_c down
% with it (3.573 -> 1.657), and Jslp/Jlcb trade nearly one-for-one, so J is flat
% across ace 0..1 and a long fit on that axis would buy ~16%. The blocking term
% is level_c, the ABR's level dependence.
%
% AND THE OBVIOUS FIX IS NOT AVAILABLE, because it is already applied:
% `tbabr_condition` sets pa.hbnl=1 UNCONDITIONALLY (tdm26.m:1536), so every ABR
% number this project has ever produced was already compressive. Setting hbnl on
% a parameter set does nothing -- measured bit-identical.
%
% What HAS never been touched is how HARD it compresses:
%     gam = cp.gm ./ (1 + hbsc * log(max(|dh|/hbmx, 1)))          micro26.m:69
%     hbsc = 0.04   hbmx = 6e-9    NEITHER IS IN parnames()
% level_c is 3.573 against a target of 5 with compression ON, so the level
% dependence is too WEAK, not missing. hbsc scales it directly; hbmx sets where
% it starts (d2mx/hbmx is currently ~21, so the knee is well below the operating
% point and lowering hbmx compresses over a wider range).
%
% tbabr_condition overrides ONLY hbnl and ihceq, so hbsc/hbmx pass through from
% pa. Verified by reading the function, but the anchor row below re-checks it:
% the first row is the untouched baseline and must reproduce 3.573.
%
% Swept at ace=+0.40, the best point on the full-objective ace axis, because that
% is where a joint fit would actually sit. maperr is a CONTROL: it comes from
% fdm26, a linear frequency-domain solve that never sees hbsc, so it must not
% move. If it does, something is wrong with the run.
t0=tic;
L=load('sweep_nch1b.mat'); pa0=L.R.pa; pa0.ace=0.40;
fprintf('\n== compression strength/threshold vs level_c, at ace=+0.40 ==\n');
fprintf('   baseline hbsc %.3g hbmx %.2g -> level_c 2.818, slope 0.5067, maperr 90.52\n', ...
        pa0.hbsc, pa0.hbmx);
fprintf('   want level_c UP toward the [3.5 6.5] band (target 5) at unchanged slope.\n\n');
fprintf('    hbsc      hbmx     slope   level_c   Jslp    Jlcb    Jsum  |  maperr(control)\n');
fprintf('   --------------------------------------------------------------------------------\n');
cfg = { pa0.hbsc, pa0.hbmx;      % anchor: must reproduce 2.818
        0.08,     pa0.hbmx;
        0.16,     pa0.hbmx;
        0.32,     pa0.hbmx;
        pa0.hbsc, 6e-10;
        pa0.hbsc, 6e-11;
        0.16,     6e-10 };
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
save('hbsweep.mat','R','cfg');
fprintf('\n   best Jslp+Jlcb: row %d, %.4f  (ace=0.40 baseline 0.1619; fitted optimum 0.1928)\n', bj, best);
fprintf('   %.1f min.  maperr must be ~90.52 on EVERY row -- it is the control.\n', toc(t0)/60);
disp('HBSWEEP_DONE');
