% PROFNEURAL -- function-level profile of one time-march condition.
%
% THE QUESTION (SN, 2026-08-02): is the neural-rate portion worth optimizing?
%
% What is already known is only STAGE level: abr_metric is 51-62% of every
% parfit26 evaluation and the only stage with an order of magnitude in it
% (nch=4: >2000 s vs score26's 681 s). What has NEVER been measured is how that
% march time splits between the MECHANICS and the NEURAL path. Without that,
% optimizing the neural code is a guess -- Amdahl decides whether it is worth
% anything at all, and if the neural share is 15% then halving it buys 7%.
%
% Reading the code first, so the profile has a prediction to falsify:
%   - hc_step (tdm26.m:362) is the neural path, called ONCE per step.
%   - Its expensive half -- the synapse -- is already SUBSAMPLED 10x:
%     `if ((nf.nihc>1) && mod(cur.stp,nf.nihc)), return; end`  (tdm26.m:376)
%     so IHC_SUBSTEPS=10 means the synapse updates every 10th step, not that the
%     IHC is oversampled. What runs EVERY step is the MET transduction: two
%     exp() over n=1401 plus a few vector ops.
%   - accel (mechanics + fluid solve) is called TWICE per step -- once directly
%     and once inside the nimp corrector loop (tdm26.m:145,152).
% PREDICTION: accel dominates and hc_step is the minority. If the profile says
% otherwise, the prediction was wrong and the neural path IS the target.
%
% ONE CONDITION, NOT THE 16-CELL GRID. tdm26('wnr1',...) runs a single
% tbabr_condition (tdm26.m:1514), the same function mix as the full grid at 1/16
% the cost. Both nch=1 and nch=3, because the split should MOVE with chamber
% count: more chambers is more fluid work, so the neural share should shrink. If
% it does not, that itself is informative.
fprintf('\n== function-level profile of one tbabr condition ==\n');
fprintf('   prediction: accel (mechanics+fluid, 2 calls/step) >> hc_step (neural, 1/step)\n');
for nch = [1 3]
    pa = modpar26(nch);
    pr = struct('fr',2,'lv',60,'pa',pa);
    profile('off'); profile('clear'); profile('on');
    t0=tic; try, evalc('S=tdm26(''wnr1'',pr,0,0);'); catch e
        fprintf('\n  nch=%d FAILED: %s\n', nch, e.message(1:min(80,end))); profile('off'); continue;
    end
    wall=toc(t0); profile('off');
    P = profile('info'); F = P.FunctionTable;
    [~,ord] = sort([F.TotalTime],'descend');
    fprintf('\n  ---- nch=%d  (%.1f s wall, profiler overhead included) ----\n', nch, wall);
    fprintf('    %-34s %8s %8s %10s\n','function','sec','%%','calls');
    tot = max(sum([F.TotalTime]),eps);
    for i = ord(1:min(14,numel(ord)))
        f = F(i);
        fprintf('    %-34s %8.2f %7.1f%% %10d\n', f.FunctionName, f.TotalTime, ...
                100*f.TotalTime/tot, f.NumCalls);
    end
    % explicit roll-up of the two candidates, by name match
    nm = {F.FunctionName};
    gi = @(s) sum([F(contains(nm,s)).TotalTime]);
    fprintf('    %s\n', repmat('-',1,64));
    fprintf('    NEURAL   hc_step   %8.2f s  %5.1f%%\n', gi('hc_step'), 100*gi('hc_step')/tot);
    fprintf('    MECH     accel     %8.2f s  %5.1f%%\n', gi('accel'),   100*gi('accel')/tot);
    fprintf('             micro26   %8.2f s  %5.1f%%\n', gi('micro26'), 100*gi('micro26')/tot);
    fprintf('             macro26   %8.2f s  %5.1f%%\n', gi('macro26'), 100*gi('macro26')/tot);
end
profile('off');
fprintf('\n  READ: the NEURAL row is the ceiling on what neural optimization can buy.\n');
fprintf('  Regardless of the split, tbabr_protocol''s 4x4 grid (tdm26.m:906-910) is 16\n');
fprintf('  INDEPENDENT conditions writing disjoint cells -- a parfor there is a pure\n');
fprintf('  scheduling change with BIT-IDENTICAL results, and cannot break a guard.\n');
disp('PROFNEURAL_DONE');
