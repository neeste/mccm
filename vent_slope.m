% DOES clvent MOVE THE FORWARD-LATENCY SLOPE?
%
% THE PROBLEM. m4_levelc found the nested m=4's one hard weakness: slope 0.126
% (no resonance) / 0.157 (clvoct 4) against a 0.413 target, where the fitted
% m=3b scaffolding sits at 0.400. A factor of 2.6-3.3 too low, on a metric that
% is a primary target of the whole TBABR+TBOAE effort. m=4's level_c turned out
% FINE (14.6 vs the scaffolding's 17.08), so the slope is the real gap.
%
% WHY clvent IS THE CANDIDATE. The two constructions BRACKET the target:
%     appended  0.597   too high
%     target    0.413
%     nested    0.157   too low
% The rebuild carried the slope past the target rather than away from it, so an
% intermediate configuration may land on it. clvent sets how strongly CL couples
% to SS and is the parameter that most directly controls how nested the model
% effectively is. vent_ss swept it for amplifier and stability and never looked
% at slope.
%
% WHAT WOULD SETTLE IT
%   slope RESPONDS to clvent -> it is a fittable lever, not a structural defect,
%                               and the fit has somewhere to go
%   slope FLAT across clvent -> the nested topology sets the slope by
%                               construction, no parameter recovers it, and m=4
%                               cannot meet the latency target in this form
%
% clvent = 3 IS A DELIBERATE REPEAT. m4_levelc measured that configuration at
% slope 0.157. If this run disagrees, the harness is wrong and nothing else in
% the table can be read -- check that row before the trend.
%
% COST each m=4 abr_metric is ~30 min, so this is ~2.5 h before contention.
% n=1401 throughout: eval_cost has not yet shown whether n=701 is FAITHFUL for
% the nested build, and a faster number that differs is a different model.
%
% References: m=3b 0.400 / maperr 499.3. m=4 default 0.157 / 525.2 / amp +56.59.

VENT = [0.5 1 3 10 30];

fprintf('\n  target slope 0.413 | m=3b scaffolding 0.400 | m=4 appended 0.597\n');
fprintf('  m=4 nested clvent=3: slope 0.157 (m4_levelc) -- REPEATED here as a check\n\n');
fprintf('  clvent | slope | level_c | maperr | amp d1  | maxRe\n');
fprintf('%s\n', repmat('-',1,60));
R = struct('v',{},'slope',{},'lc',{},'maperr',{},'amp',{});
for v = VENT
    pa = modpar26(4); pa.clvent = v;
    sl = NaN; lc = NaN; mp = NaN; ag = NaN; mr = NaN;
    try
        evalc('m = abr_metric(pa,false);'); sl = m.slope; lc = m.level_c;
    catch e
        fprintf('  %6.1f | abr_metric FAILED: %s\n', v, e.message);
    end
    try
        S = score26(pa,'fast',false); mp = S.maperr; ag = S.amp_gain; mr = S.maxRe;
    catch
    end
    fprintf('  %6.1f | %5.3f | %7.2f | %6.1f | %+6.2f | %+7.1f\n', v, sl, lc, mp, ag, mr);
    R(end+1).v=v; R(end).slope=sl; R(end).lc=lc; R(end).maperr=mp; R(end).amp=ag; %#ok<SAGROW>
    save('vent_slope.mat','R');
end
fprintf(['\n  CHECK clvent=3 AGAINST 0.157 FIRST. Then read the trend: if slope\n' ...
         '  climbs toward 0.413 anywhere in this range, clvent is a real lever and\n' ...
         '  m=4''s latency gap is a tuning problem. If it sits near 0.15 throughout,\n' ...
         '  the nested topology fixes the slope by construction and no amount of\n' ...
         '  fitting this parameter will reach the target.\n' ...
         '  Watch amp d1 alongside: a clvent that fixes the slope but pushes the\n' ...
         '  amplifier out of 40-60 has moved the problem rather than solved it.\n']);
disp('VENT_SLOPE_DONE');
