% GPSWEEP -- fine gampro-grade sweep from the NEW nch=1 best fit (sweep_nch1b).
%
% WHY. The nch=1 rerun finished at J=0.1928 with map/gain/osc terms all EXACTLY
% zero, so J is now 100% the slope term: |0.606 - 0.413|. 150 evaluations of 15
% impedance/chsz parameters moved the slope only 0.663 -> 0.606. The fit vector
% is exhausted on the one quantity that still costs anything.
%
% pa.gampro (place-dependent CA gain) is ones(n,1) in every parameter set and
% has NEVER been fitted -- it is not in parnames(), so it is not in the fit
% vector at all. The 2026-07-22 lever sweep measured enormous slope leverage at
% EXACTLY constant maperr (150.5 at every g), which is the one property the 15
% fitted parameters do not have.
%
% THE CATCH, and the reason this is measured before anything is fitted: the
% direction we need is the direction that went unstable. That sweep, from a
% DIFFERENT baseline, read osc -15 at g=+0.4 and osc +16659 at g=-0.4, and we
% need g<0 to LOWER the slope. The present baseline sits at osc -56.6, i.e.
% only 16.6 dB of headroom before parfit26's own osc guard (maxRe_osc > -40)
% starts charging. So the question is not "does gampro move the slope" -- that
% is known -- but "does it reach 0.413 before the osc guard bites".
%
% Measured exactly as jointobj measures them (abr_metric for slope, score26
% 'fast' for maperr/maxRe/maxRe_osc/amp_gain), so the columns are directly
% comparable to a parfit26 run with wgain=0.01.
%
% NOTE for whoever fits this next: jointobj builds pa from base=modpar26(nch)
% plus the pv vector, so a gampro carried on a warm start would be SILENTLY
% RESET to ones(n,1) on every evaluation. Fitting it requires adding it to the
% parameter vector (or pinning it, which holds it out of the fit). Do not
% assume a warm start carries it.
t0=tic;
L=load('sweep_nch1b.mat'); pa0=L.R.pa;
n=pa0.n; xf=((0:n-1)')/(n-1);
% GRID REFINED after the first pass (coarse grid 0 / -0.05 / ... / -0.30 killed
% once the second point answered the question). Measured leverage is ~6.7 slope
% per unit g -- SIX TIMES the 2026-07-22 estimate of ~1.05 -- so g=-0.05 already
% overshoots the 0.413 target down to 0.272, while maxRe_osc goes -56.6 -> +566.8.
% maperr held at EXACTLY 90.97, so the CF-map-free property is confirmed and the
% only live question is where the stability wall sits. Everything below -0.05 is
% known-unstable and carries no information. Sample where the target is crossed.
% Anchor row g=0 kept deliberately: a decomposition without a row whose answer is
% already known is how a wrong Df/Dq conclusion got through on 2026-07-29.
gg=[0 -0.010 -0.015 -0.020 -0.025 -0.030 -0.035 -0.040];
fprintf('\n== gampro grade sweep, gampro=exp(g*(xf-0.5)), from sweep_nch1b (nch=1) ==\n');
fprintf('   baseline g=0: slope 0.606  maperr 90.97  amp +54.54  osc -56.6\n\n');
fprintf('     g     slope   maperr   amp_gain   maxRe   maxRe_osc |  Jslope   Jmap    Jgain   Josc   =  J\n');
fprintf('   -----------------------------------------------------------------------------------------------\n');
R=struct('g',num2cell(gg));
for j=1:numel(gg)
    g=gg(j); pa=pa0; pa.gampro=exp(g*(xf-0.5));
    sl=NaN; mp=NaN; ag=NaN; mr=NaN; mo=NaN;
    try
        m=abr_metric(pa,false);
        if (m.ok), sl=m.slope; end
    catch e
        fprintf('   %+5.2f  abr_metric FAILED: %s\n', g, e.message(1:min(60,end)));
    end
    try
        S=score26(pa,'fast',false);
        mp=S.maperr; mr=S.maxRe; mo=S.maxRe_osc;
        if (isfield(S,'amp_gain')&&~isempty(S.amp_gain)), ag=double(S.amp_gain(1)); end
    catch e
        fprintf('   %+5.2f  score26 FAILED: %s\n', g, e.message(1:min(60,end)));
    end
    Js=abs(sl-0.413); Jp=0.001*max(0,mp-105); Jg=0.01*max(0,40-ag);
    Jo=0.005*max(0,mo+40); J=Js+Jp+Jg+Jo;
    fprintf('   %+5.2f  %6.3f  %7.2f   %+7.2f  %7.1f    %7.1f  | %6.4f  %6.4f  %6.4f  %6.4f  = %6.4f\n', ...
            g, sl, mp, ag, mr, mo, Js, Jp, Jg, Jo, J);
    R(j).slope=sl; R(j).maperr=mp; R(j).amp=ag; R(j).maxRe=mr; R(j).osc=mo; R(j).J=J;
end
save('gpsweep.mat','R','gg');
fprintf('\n   %.1f min.  Read: does slope reach 0.413 while maxRe_osc stays below -40?\n', toc(t0)/60);
disp('GPSWEEP_DONE');
