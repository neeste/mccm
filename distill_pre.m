function distill_pre
%DISTILL_PRE  Pre-flight for the nch=1,dof=3 distillation fit.
%
% The design note: "The 1-chamber model is a distillation of the 3-D physics,
% obtained by fitting an effective partition impedance to reproduce the
% 3-chamber model's response, NOT by hand-tuning against animal data."
%
% TARGET is the 3-chamber model's cfmap -- peak place per frequency for the BM
% and the hair bundle, which is exactly what the note names. cfmap_analysis is
% reused rather than reinvented: it exists because fdmod23/tuning_err are local
% to fdm26 and four external harnesses failed trying to reach them.
%
% DECISION 2 IS UNRESOLVED and governs the target: linear FD response only, or
% also time-domain level dependence? This builds the LINEAR FD target, which is
% a necessary subset of either answer -- if the 1-chamber form cannot match the
% linear response, the harder target is moot.
%
% TWO PRE-FLIGHT RULES, both learned expensively in this project:
%   1. Evaluate the objective at the START POINT before spending a fit. A
%      7-hour fit was lost when the metric rejected its own starting model.
%   2. Check it RESPONDS, not just that it evaluates. The SS=0.10 refit ran 153
%      evaluations against J=0.0000 because maptol hinged the map term out of
%      existence -- a constant objective evaluates perfectly well.

FL = 500*2.^(-1:5);
fprintf('\n  TARGET: 3-chamber cfmap at %d frequencies\n', numel(FL));
t0=tic; T = fdm26(struct('cfmap',1,'pa',modpar26(3),'flst',FL)); tt=toc(t0);
fprintf('      xpk_bm: %s\n', num2str(T.xpk_bm,'%.4f '));
fprintf('      xpk_hb: %s\n', num2str(T.xpk_hb,'%.4f '));
fprintf('      (%.1f s per target evaluation)\n', tt);

p0 = modpar26(1); p0.dof = 3;
fprintf('\n  START: nch=1, dof=3\n');
t0=tic; J0 = distill_obj(p0, T, FL); te=toc(t0);
fprintf('      J = %.6g   (%.1f s per objective evaluation)\n', J0, te);
if (~isfinite(J0))
    fprintf('      *** NON-FINITE at the start point -- the fit would begin from\n');
    fprintf('          nowhere. This is the 7-hour failure. STOP. ***\n'); return
end

% --- does it RESPOND? perturb each candidate parameter -------------------
flds = {'k5o','k5e','r5o','m5o','m5e','k1o','k1e','m1o'};
fprintf('\n  RESPONSIVENESS: does J move when each parameter moves?\n');
fprintf('      param |    J(+2%%)    |   dJ/J0    | live\n');
fprintf('      %s\n', repmat('-',1,46));
nlive = 0;
for i = 1:numel(flds)
    f = flds{i}; p = p0;
    if (~isfield(p,f)), fprintf('      %-5s | (absent)\n', f); continue; end
    p.(f) = p.(f) * 1.02;
    J = distill_obj(p, T, FL);
    d = (J - J0)/max(abs(J0),eps);
    live = abs(d) > 1e-9;
    nlive = nlive + live;
    lv = 'no '; if (live), lv = 'YES'; end
    fprintf('      %-5s | %12.6g | %+10.3e | %s\n', f, J, d, lv);
end
fprintf('\n      %d of %d parameters move the objective\n', nlive, numel(flds));
if (nlive == 0)
    fprintf('      *** OBJECTIVE IS CONSTANT -- nothing to fit. STOP. ***\n');
else
    fprintf('      GO: objective is finite and responsive.\n');
end
fprintf(['\n  NOTE this fits the 1-chamber model to the 3-CHAMBER MODEL, not to\n' ...
         '  data. A good J means faithful distillation, not a good cochlea --\n' ...
         '  the mechanism model carries whatever error it has into the target.\n']);
end

function J = distill_obj(pa, T, FL)
% Mismatch between this model's cfmap and the 3-chamber target. Place error in
% OCTAVES-equivalent terms is not available, so normalized place is used
% directly; both BM and hair-bundle peaks are weighted equally, per the design
% note naming both.
try
    R = fdm26(struct('cfmap',1,'pa',pa,'flst',FL));
catch
    J = Inf; return
end
eb = R.xpk_bm - T.xpk_bm;
eh = R.xpk_hb - T.xpk_hb;
ok = isfinite(eb) & isfinite(eh);
if (~any(ok)), J = Inf; return; end
J = sqrt(mean(eb(ok).^2) + mean(eh(ok).^2));
end
