% DOES PROXIMITY TO CRITICALITY DRIVE THE NON-CONVERGENCE?
%
% SN's second lead, and it unifies with the first. Near-critical operation means
% the active force nearly cancels the passive damping, so effective Q is high
% and the system is STIFF. A stiff system is exactly where ONE corrector
% iteration (pa.nimp=1) fails to converge. So low net damping and low nimp are
% not separate problems: the implicit solve's convergence RATE degrades as the
% model approaches criticality.
%
% EXISTING EVIDENCE, from conv_slope read the right way:
%     model  maxRe   slope across n           level_c across n
%     m=3    +19.3   0.117 / 0.400 / 0.255    31.96 / 17.08 / 29.92
%     m=4     +2.3   0.333 / 0.157 / 0.140    13.97 / 14.60 / 15.11
% m=3b sits 8x closer to criticality and is far worse conditioned. I first read
% that as a surprising reversal; under this hypothesis it is expected.
%
% THIS ALSO QUALIFIES A CLAIM I HAVE LEANED ON. I have cited m=4's low maxRe as
% evidence of a well-behaved amplifier. The same margin that makes it stable
% makes it numerically tractable, so part of "m=4 is better conditioned" may
% just be "m=4 is further from the operating point the cochlea uses".
%
% THE TEST holds the model fixed and varies ONLY gam (the NDR multiplier, which
% scales the active force and therefore the net damping), measuring grid
% sensitivity at each setting. gam=1 is native; lower gam means more net damping
% and greater distance from criticality.
%
% SENSITIVITY is reported as the SPREAD of each metric across n=701/1401/2801,
% which is the quantity in question -- not the metric values themselves.
%   spread SHRINKS as gam falls -> criticality drives the non-convergence, and
%     the model must be integrated differently (higher nimp, smaller dt) to be
%     trusted anywhere near its real operating point
%   spread FLAT in gam           -> stiffness is not the mechanism, and the
%     dt/dx confound stays the leading explanation
%
% m=3 is used because abr_metric is ~102 s there against ~2290 s at m=4, and
% because m=3 showed the worst grid behaviour. maxRe is reported so the actual
% distance from criticality is visible rather than assumed from gam.

GAM = [1.0 0.7 0.4];
NN  = [701 1401 2801];

fprintf('\n  Varying ONLY gam (net damping). Spread across n is the quantity of\n');
fprintf('  interest, not the values. At gam=1 the spreads were:\n');
fprintf('    slope 0.283 (0.117-0.400) | level_c 14.88 | maperr 388.6\n\n');
fprintf('  gam  | maxRe   | slope: 701/1401/2801        | spread | maperr spread\n');
fprintf('%s\n', repmat('-',1,76));
for g = GAM
    sl = nan(1,numel(NN)); mp = nan(1,numel(NN)); mr = NaN;
    for i = 1:numel(NN)
        pa = modpar26(3); pa.gam = g;
        if (NN(i) ~= pa.n)
            try, pa = setn(pa, NN(i)); catch, continue; end
        end
        try, evalc('m = abr_metric(pa,false);'); sl(i) = m.slope; catch, end
        try
            S = score26(pa,'fast',false); mp(i) = S.maperr;
            if (NN(i)==1401), mr = S.maxRe; end
        catch
        end
    end
    ssl = NaN; smp = NaN;
    if (sum(isfinite(sl))>=2), ssl = max(sl(isfinite(sl)))-min(sl(isfinite(sl))); end
    if (sum(isfinite(mp))>=2), smp = max(mp(isfinite(mp)))-min(mp(isfinite(mp))); end
    fprintf('  %.1f  | %+7.1f | %6.3f %6.3f %6.3f | %6.3f | %8.1f\n', ...
        g, mr, sl(1), sl(2), sl(3), ssl, smp);
end
fprintf(['\n  If the spreads shrink as gam falls, criticality drives the\n' ...
         '  non-convergence. That would mean the model can only be trusted at\n' ...
         '  operating points well away from the one it is meant to represent,\n' ...
         '  unless the integration is changed -- and it would make nimp and dt\n' ...
         '  prerequisites for fitting rather than refinements.\n' ...
         '  Note gam also changes the MODEL, so the metric values are not\n' ...
         '  comparable across rows; only the SPREADS are.\n']);
disp('DAMP_CONV_DONE');
