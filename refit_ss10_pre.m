% PRE-FLIGHT FOR THE SS=0.10 REFIT. Evaluate the objective at the START POINT
% before spending a fit.
%
% THIS IS A HARD RULE IN THIS PROJECT and it is written in blood: a 7-hour fit
% was lost when tiptail_metric's 0.6-16 kHz window REJECTED ITS OWN STARTING
% MODEL, so the search began from an infeasible point and never recovered. One
% cheap evaluation up front would have caught it.
%
% WHAT IS BEING SET UP
%   pin      chsz = [0.95 0.10 1.00]. opts.pin forces the field onto BOTH base
%            and warm start so setpar_l cannot overwrite it mid-fit.
%   fitidx   DEFAULT MINUS 31,32,33. Those three indices ARE the chamber sizes;
%            leaving them in would refit chsz and undo the whole point of the
%            exercise. This is the single most important line here.
%   warm     modpar26(3) with the pinned chsz, NOT the default
%            refit_c3_map.mat. That file exists and would otherwise be used, but
%            it carries different parameters from modpar26c3's inline values --
%            and the inline values are what produced the 406.8 measurement this
%            refit is trying to improve on. Starting elsewhere would make the
%            before/after incomparable.
%
% BASELINE TO BEAT: m=3b at [0.95 0.10 1.00] scores maperr 406.8 UNFITTED,
% against the fitted stock [0.95 0.05 1.00] at 499.3. The question is how much
% of the remaining error is the area mismatch and how much is genuine.
%
% maxfe=1 so this returns after the start evaluation. It also times one
% evaluation, which is what sizes the real budget.

pa = modpar26(3);
pa.chsz = [0.95 0.10 1.00];

opts = struct();
opts.pin    = struct('chsz', [0.95 0.10 1.00]);
opts.fitidx = [1 2 3 4 5 7 8 9 10 11 13 20 21];   % default MINUS 31,32,33 (chsz)
opts.warm   = pa;
opts.maxfe  = 1;
opts.out    = 'refit_ss10_pre.mat';

fprintf('\n  fitidx = [%s]\n', strtrim(sprintf('%d ', opts.fitidx)));
fprintf('  chsz pinned to [%s] (sum %.2f)\n', ...
    strtrim(sprintf('%.2f ', opts.pin.chsz)), sum(opts.pin.chsz));
fprintf('  UNFITTED baseline at these areas: maperr 406.8\n');
fprintf('  fitted stock [0.95 0.05 1.00]:    maperr 499.3\n\n');

t0 = tic;
R = parfit26(3, opts);
el = toc(t0);

fprintf('\n  --- PRE-FLIGHT SUMMARY ---\n');
fprintf('  elapsed %.1f s for a ~1-evaluation run\n', el);
if (isfield(R,'Rf') && isfield(R.Rf,'maperr'))
    fprintf('  maperr at start : %.1f\n', R.Rf.maperr);
end
if (isfield(R,'pa'))
    fprintf('  chsz after      : [%s]  <== MUST still be 0.95 0.10 1.00\n', ...
        strtrim(sprintf('%.3f ', R.pa.chsz)));
end
fprintf(['\n  GO if the start objective is finite, maperr is near 406.8, and chsz\n' ...
         '  survived unchanged. NO-GO if chsz moved (the pin or fitidx is wrong),\n' ...
         '  or if the objective is Inf/NaN (the start point is infeasible and the\n' ...
         '  search would begin from nowhere -- the exact 7-hour failure).\n']);
disp('REFIT_SS10_PRE_DONE');
