% CAN m=3b REPLACE m=3 WITHOUT LOSS OF GENERALITY?
% The click identity already passed exactly (d1/d2/ped bit-identical). Two
% things remain:
%  (1) coupeig reports a LARGER EIGENVALUE SET for m=3b (3 DOFs, so d3's modes
%      are included). Those modes are decoupled and cannot affect d1/d2, but
%      maxRe is a max over ALL modes, so it could differ for a non-physical
%      reason. If it matches, consolidation is clean on stability too.
%  (2) COST. State goes 2n -> 3n and coupeig's dense eig scales ~N^3, so the
%      stability check should cost ~3.4x more.
% fdm26/maperr is unaffected either way: it has no d3 for a 3-chamber model.
cfg = { 'm=3 ', modpar26(3); 'm=3b', modpar26c3b };
fprintf('\n  model | maperr | maxRe     | maxRe_osc | range mono | contr | amp d1 | sec\n');
fprintf('%s\n', repmat('-',1,80));
S = cell(1,2);
for i = 1:2
    t0 = tic; S{i} = score26(cfg{i,2}, 'fast', false); sec = toc(t0);
    fprintf('  %-5s | %6.1f | %+9.1f | %+9.1f | %5.2f %-4s | %5.1f | %+6.2f | %4.0f\n', ...
        cfg{i,1}, S{i}.maperr, S{i}.maxRe, S{i}.maxRe_osc, S{i}.bf_range, ...
        S{i}.bf_mono, S{i}.contrast, S{i}.amp_gain, sec);
end
a=S{1}; b=S{2};
fprintf('\n  differences (m=3b minus m=3):\n');
fprintf('    maperr    %+.3f\n', b.maperr-a.maperr);
fprintf('    maxRe     %+.3f\n', b.maxRe-a.maxRe);
fprintf('    maxRe_osc %+.3f\n', b.maxRe_osc-a.maxRe_osc);
fprintf('    amp d1    %+.3f dB\n', b.amp_gain-a.amp_gain);
fprintf('    contrast  %+.3f dB\n', b.contrast-a.contrast);
if (abs(b.maxRe-a.maxRe)<1e-6 && abs(b.maperr-a.maperr)<1e-6 && abs(b.amp_gain-a.amp_gain)<1e-6)
    fprintf('\n  CLEAN: every scored quantity matches. m=3b can replace m=3.\n');
else
    fprintf('\n  DIFFERS: consolidation would change reported values, see above.\n');
end
disp('CONSOL_CHECK_DONE');
