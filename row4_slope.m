% SLOPE FOR THE BEST m=4 AREAS -- the metric that has caught this twice.
%
% ROW 4 of dual_bisect, chsz [1.0 0.10 0.85 0.05], is the best m=4 found:
%   maperr 466.3  amp +40.48  maxRe +3.3  chi 8.9
% It beats the FITTED m=3b scaffolding (499.3) by 6.6% while itself being
% UNFITTED, with the amplifier inside SN's 40-60 band. On score26 fast metrics
% it is the strongest configuration in the project.
%
% WHY THAT IS NOT ENOUGH. score26 fast does not compute slope or level_c. The
% SS=0.10 result at m=3 looked equally good on those same metrics -- maperr
% 499.3 -> 406.8, the only clean fold in the project, contrast 9.2 -> 14.9 --
% until the forward-latency slope came in at 0.185 against a 0.413 target,
% having been 0.400 at SS=0.05. SS area trades CF map against slope, and row 4
% uses exactly that lever. Its map win may be bought the same way.
%
% THE QUESTION: does m=4 pay m=3's price for SS=0.10?
%   slope near 0.400  -> the nested construction does NOT pay it, and row 4 is
%                        a genuine advance rather than a lever being pulled
%   slope near 0.185  -> same trade as m=3, and row 4's map win is not free
%
% ROW 2 is the other map-improving variant from dual_bisect: CL doubled instead
% of SS (maperr 495.1, range 5.87, fold 0.43, but amp +70.76 is above band).
% Including it separates the two carve levers on the slope axis -- if CL
% doubling is slope-neutral where SS doubling is not, that is the more useful
% lever even though its amplifier needs retuning.
%
% REFERENCES  m=3b stock (fitted): slope 0.400, level_c 17.08, maperr 499.3
%             m=3  at SS=0.10    : slope 0.185, level_c 36.75, maperr 406.8
%             target slope 0.413; level_c target 5 but NOTHING reaches its band,
%             the fitted scaffolding included, so read level_c against 17.08.
%
% m=4's slope for the CURRENT areas comes from the m4_levelc run, so it is not
% repeated here.

cfg = { 'row4 SS=0.10 [1.0 .10 .85 .05]', [1.00 0.10 0.85 0.05]
        'row3 CL=0.10 [1.0 .05 .85 .10]', [1.00 0.05 0.85 0.10] };

fprintf('\n  m=3b stock : slope 0.400  level_c 17.08  maperr 499.3\n');
fprintf('  m=3 SS=0.10: slope 0.185  level_c 36.75  maperr 406.8  <- the price\n');
fprintf('  target slope 0.413\n\n');
fprintf('  config                          | slope | level_c | %%/dB  | maperr | amp d1\n');
fprintf('%s\n', repmat('-',1,78));
R = struct('nm',{},'slope',{},'lc',{},'maperr',{});
for i = 1:size(cfg,1)
    cz = cfg{i,2};
    if (abs(sum(cz)-2) > 1e-12)
        fprintf('  %-31s | SUM %.4f ~= 2 -- SKIPPED\n', cfg{i,1}, sum(cz)); continue
    end
    pa = modpar26(4); pa.chsz = cz;
    sl = NaN; lc = NaN; pdb = NaN; mp = NaN; ag = NaN;
    try
        evalc('m = abr_metric(pa, false);');
        sl = m.slope; lc = m.level_c;
        if (isfinite(lc) && lc > 0), pdb = 100*(lc^(1/100)-1); end
    catch e
        fprintf('  %-31s | abr_metric FAILED: %s\n', cfg{i,1}, e.message);
    end
    try
        S = score26(pa, 'fast', false); mp = S.maperr; ag = S.amp_gain;
    catch
    end
    fprintf('  %-31s | %5.3f | %7.2f | %5.2f | %6.1f | %+6.2f\n', ...
        cfg{i,1}, sl, lc, pdb, mp, ag);
    R(end+1).nm=cfg{i,1}; R(end).slope=sl; R(end).lc=lc; R(end).maperr=mp; %#ok<SAGROW>
    save('row4_slope.mat','R');
end
fprintf(['\n  Slope near 0.400 means m=4 does NOT pay m=3''s price and row 4 is a\n' ...
         '  real advance. Near 0.185 means the same trade, and row 4''s maperr\n' ...
         '  466.3 is bought with the forward latency -- in which case the honest\n' ...
         '  comparison against m=3b is at matched slope, not matched areas.\n']);
disp('ROW4_SLOPE_DONE');
