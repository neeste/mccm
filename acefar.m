% ACEFAR -- push ace past +0.73 and find where the OBJECTIVE bottoms out.
%
% acecheck established ace as a real, sign-favorable lever: slope 0.606 -> 0.474
% across ace -0.399 -> +0.73, with maxRe_osc IMPROVING 34.7 dB and maperr dipping
% to 79.17 (better than the 90.97 project best) before rising again. Nothing else
% in the 16-parameter atlas does that; every other parameter pays for slope with
% stability.
%
% WHY THIS MATTERS MORE THAN IT LOOKS. ace is pv index 20 and was ALREADY in the
% fitted set. The 150-evaluation run moved it 0.20% -- from -0.400000 to
% -0.399204 -- while the sweep says it wants to move ~440%. So the fit was never
% short of levers, it was short of SEARCH: 15-dimensional Nelder-Mead spends 16
% evaluations just building its simplex, and the run's J flatlined at 0.1928 from
% iteration 74 to 84, i.e. the simplex collapsed before this coordinate was ever
% explored. Scored against the real objective, ace=+0.73 alone is J ~ 0.061
% versus the fitted 0.1928.
%
% WHAT STOPS IT is NOT stability, which keeps improving. It is maperr climbing
% back toward maptol=105 and amp_gain falling toward gainmin=40 (100.20 and
% +41.96 at ace=0.73, both close). Those are exactly the constraints the other 14
% parameters exist to repair, which is the argument for a JOINT fit warm-started
% from here rather than a hand-set value.
%
% J is scored exactly as jointobj scores it, so the minimum found here is
% directly comparable to the 0.1928 on record.
t0=tic;
L=load('sweep_nch1b.mat'); pa0=L.R.pa;
vv=[0.73 1.00 1.35 1.70 2.10];
fprintf('\n== ace continued, from sweep_nch1b (nch=1) ==\n');
fprintf('   on record: ace -0.3992 -> J 0.1928 (slope 0.60581, the fitted optimum)\n');
fprintf('              ace +0.7300 -> J 0.0610 (slope 0.47398), never visited by the fit\n\n');
fprintf('     ace      slope     maperr    amp     osc    |  Jslope   Jmap    Jgain   Josc  =    J\n');
fprintf('   ----------------------------------------------------------------------------------------\n');
best=Inf; bace=NaN; R=struct('ace',num2cell(vv));
for j=1:numel(vv)
    pa=pa0; pa.ace=vv(j);
    sl=NaN; mp=NaN; ag=NaN; mo=NaN;
    try, m=abr_metric(pa,false); if (m.ok), sl=m.slope; end, catch, end
    try
        S=score26(pa,'fast',false); mp=S.maperr; mo=S.maxRe_osc;
        if (isfield(S,'amp_gain')&&~isempty(S.amp_gain)), ag=double(S.amp_gain(1)); end
    catch, end
    Js=abs(sl-0.413); Jp=0.001*max(0,mp-105); Jg=0.01*max(0,40-ag);
    Jo=0.005*max(0,mo+40); J=Js+Jp+Jg+Jo;
    if (isfinite(J) && J<best), best=J; bace=vv(j); end
    fprintf('   %+7.4f  %8.5f  %8.2f  %+6.2f  %7.1f  | %6.4f  %6.4f  %6.4f  %6.4f = %6.4f\n', ...
            vv(j), sl, mp, ag, mo, Js, Jp, Jg, Jo, J);
    R(j).slope=sl; R(j).maperr=mp; R(j).amp=ag; R(j).osc=mo; R(j).J=J;
end
save('acefar.mat','R','vv','best','bace');
fprintf('\n   best on this axis: ace %+.4f  J %.4f   (fitted optimum was J 0.1928)\n', bace, best);
fprintf('   %.1f min.  This point is the warm start for the long joint fit.\n', toc(t0)/60);
disp('ACEFAR_DONE');
