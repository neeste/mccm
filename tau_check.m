% IS dt SMALL COMPARED TO THE RESONATOR TIME CONSTANTS?
%
% SN's point, and it is the quantitative statement the last two runs were
% circling: r/k has units of TIME and sets how fast a displacement can change.
% The corrector in tdm26's predictor-corrector step converges at a rate governed
% by dt/tau, so if dt is a large fraction of tau, ONE correction (pa.nimp=1) is
% nowhere near enough -- and no amount of SPATIAL refinement fixes it.
%
% THAT LAST POINT INVALIDATES MY OWN CONVERGENCE TEST. conv_slope refined n at
% FIXED dt=2e-6 and found both models unconverged. I read that as a dt/dx
% confound. But if the limiting error is temporal -- set by dt/tau -- then
% refining n could never have converged anything, and the test was asking a
% question that had no answer. Six hours of runs to learn that.
%
% HOW THE THREE LEADS UNIFY
%   tau = r/k          the direct measure
%   nimp               the corrector needs more iterations when dt/tau is large
%   net damping        near-critical operation cancels r, so tau -> 0 and
%                      dt/tau -> infinity. The amplifier's own mechanism is what
%                      makes the integration hard, which is why this survived
%                      unnoticed: at low gain nimp=1 is adequate.
%
% This is pure arithmetic on cp -- no time march, no fitting. One cheap tdm26
% call to obtain cp on the model's own grid rather than guessing the x spacing
% (the apical point has k1 -> 0, so a hand-built grid would misreport it).
%
% RULE OF THUMB for a trapezoidal corrector: dt/tau below ~0.1 is comfortable,
% ~0.5 means one iteration leaves a large residual, above 1 means the step
% cannot track the relaxation at all.

for nch = [3 4]
    pa = modpar26(nch);
    fprintf('\n  ===== m=%d, dt = %.2e s, nimp = %d =====\n', nch, pa.dt, pa.nimp);
    try
        evalc('S = tdm26(0,pa,0,0);'); cp = S.cp;
    catch e
        fprintf('  could not obtain cp: %s\n', e.message); continue
    end
    dt = pa.dt;
    nm = {}; tau = {};
    nm{end+1}='d1  r1/k1';            tau{end+1} = cp.r1 ./ max(cp.k1,eps);
    nm{end+1}='d2  (r2+r3)/(k2+k3)';  tau{end+1} = (cp.r2+cp.r3) ./ max(cp.k2+cp.k3,eps);
    if (isfield(cp,'k5') && ~isempty(cp.k5))
        nm{end+1}='d3  r5/k5';        tau{end+1} = cp.r5 ./ max(cp.k5,eps);
    end
    if (isfield(cp,'clvm') && ~isempty(cp.clvm) && any(cp.clvk(:)~=0))
        nm{end+1}='vent clvr/clvk';   tau{end+1} = cp.clvr ./ max(cp.clvk,eps);
    end
    fprintf('  element              | tau min s | tau max s | dt/tau MAX | verdict\n');
    fprintf('  %s\n', repmat('-',1,68));
    for i = 1:numel(nm)
        t = tau{i}; t = t(isfinite(t) & t>0);
        if (isempty(t)), fprintf('  %-20s | (no finite tau)\n', nm{i}); continue; end
        rmax = dt/min(t);
        if (rmax < 0.1),      v = 'comfortable';
        elseif (rmax < 0.5),  v = 'MARGINAL';
        elseif (rmax < 1),    v = 'BAD';
        else,                 v = 'UNRESOLVED';
        end
        fprintf('  %-20s | %9.2e | %9.2e | %10.2f | %s\n', ...
            nm{i}, min(t), max(t), rmax, v);
    end
end
fprintf(['\n  dt/tau is worst at the BASE, where k is largest and tau smallest.\n' ...
         '  Anything at MARGINAL or worse means nimp=1 leaves a residual that\n' ...
         '  spatial refinement cannot remove -- and that the fix is more corrector\n' ...
         '  iterations or a smaller dt, not a finer grid.\n' ...
         '  NOTE these are PASSIVE time constants. Near criticality the active\n' ...
         '  force cancels part of r, so the EFFECTIVE tau is smaller and dt/tau\n' ...
         '  worse than shown -- the amplifier makes its own integration harder.\n']);
disp('TAU_CHECK_DONE');
