% HBPROBE -- does MET compression free the level_c term that is blocking ace?
%
% THE CHAIN THIS TESTS.
%   1. ace reaches the ABR slope target (0.6058 -> 0.3856, crossing 0.413).
%   2. What stops it is NOT stability -- maxRe_osc IMPROVES all the way to -107.8
%      -- but the level_c band term: Jslp and Jlcb trade nearly one-for-one along
%      the ace axis, so J is flat (0.1711/0.1619/0.1642/0.1684 across ace 0..1).
%   3. level_c IS the ABR's level dependence, i.e. a COMPRESSION observable.
%   4. pa.hbnl = 0 in every parameter set, so the amplifier gain is level
%      INDEPENDENT. micro26.m:61 is the only place compression enters.
%   5. So the term blocking ace penalizes the model for something the linear
%      configuration cannot do. Turn compression on and level_c should move.
%
% FIRST CHECK IS WHETHER IT ENGAGES AT ALL. The compression law tests |dh|
% against pa.hbmx (=6e-9): gam/(1 + hbsc*log(max(|dh|/hbmx,1))). tdm26.m:1547
% warns explicitly that if d2mx << hbmx the gain never compresses and the model
% runs LINEAR regardless of the flag, and exposes dgn.ratio = d2mx/hbmx for
% exactly this. A ratio below 1 means hbnl=1 changed NOTHING and any apparent
% level_c shift is noise. Report it first, believe nothing else until it clears.
%
% NOTE score26's maperr comes from fdm26, a frequency-domain LINEAR solve, so
% maperr should be INVARIANT to hbnl. That invariance is a free control: if
% maperr moves, hbnl is reaching something it should not, and the run is wrong.
t0=tic;
L=load('sweep_nch1b.mat'); pa0=L.R.pa;
aces=[-0.3992 0.40 1.00];
fprintf('\n== hbnl probe: does compression free level_c? ==\n');
fprintf('   band [3.5 6.5] wlcb 0.10.  hbmx %.2g, hbsc %.3g, mmeq %d\n\n', ...
        pa0.hbmx, pa0.hbsc, subsref_d(pa0,'mmeq',1));
fprintf('     ace   hbnl   slope   level_c   Jslp    Jlcb   |  maperr   d2mx/hbmx  compresses?\n');
fprintf('   ---------------------------------------------------------------------------------\n');
R=struct();  k=0;
for a=aces
    for hb=[0 1]
        pa=pa0; pa.ace=a; pa.hbnl=hb;
        sl=NaN; lc=NaN; mp=NaN; rt=NaN;
        try, m=abr_metric(pa,false); if (m.ok), sl=m.slope; lc=m.level_c; end, catch, end
        try, S=score26(pa,'fast',false); mp=S.maperr; catch, end
        try
            D=tdm26('wnr1',struct('fr',2,'lv',60,'pa',pa),0,0);
            if (isstruct(D)&&isfield(D,'dgn')&&isfield(D.dgn,'ratio')), rt=D.dgn.ratio; end
        catch, end
        Js=abs(sl-0.413); Jl=0.10*(max(0,3.5-lc)+max(0,lc-6.5));
        cw='--'; if (isfinite(rt)), if (rt>1), cw='YES'; else, cw='no (linear)'; end, end
        fprintf('   %+6.3f    %d   %7.4f  %7.3f  %6.4f  %6.4f  | %7.2f  %9.3f  %s\n', ...
                a, hb, sl, lc, Js, Jl, mp, rt, cw);
        k=k+1; R(k).ace=a; R(k).hbnl=hb; R(k).slope=sl; R(k).level_c=lc;
        R(k).maperr=mp; R(k).ratio=rt;
    end
end
save('hbprobe.mat','R','aces');
fprintf('\n   %.1f min.\n', toc(t0)/60);
fprintf('   READ: (a) d2mx/hbmx > 1 or nothing happened; (b) level_c RISES toward 3.5\n');
fprintf('         with hbnl=1; (c) maperr UNCHANGED (fdm26 is linear -- a control).\n');
disp('HBPROBE_DONE');
function v=subsref_d(s,f,d), if (isfield(s,f)), v=s.(f); else, v=d; end, end
