% level_c FOR m=4 -- the metric score26 fast never computes.
%
% WHY THIS IS BEING RUN. Every m=4 number quoted today came from score26 fast,
% which reports maperr, amp, maxRe, range, fold and contrast and NOT level_c.
% That blind spot just invalidated the SS=0.10 result: it looked like an
% improvement on six metrics (maperr 499.3 -> 406.8, the only clean fold in the
% project, contrast 9.2 -> 14.9) until the refit pre-flight reported level_c =
% 36.75 against a target of 5. There is no reason to assume m=4 is safe on an
% axis nothing has checked.
%
% level_c maps to per-dB latency shift as %/dB = 100*(level_c^(1/100)-1):
%   level_c 5     -> 1.62 %/dB   the data
%   band [3.5 6.5]-> 1.26-1.88 %/dB
%   level_c 36.75 -> 3.67 %/dB   SS=0.10, more than twice the data
%
% ROW 1 IS THE CONTROL and it matters as much as the m=4 rows. m=3b stock is the
% FITTED scaffolding, so its level_c should land in band. If it does not, then
% either the metric is not measuring what the objective's target implies or the
% fitted parameters have drifted, and every reading below is uninterpretable
% until that is resolved. Do not read rows 2-4 before checking row 1.
%
% ROWS 2-4 separate the construction from the resonance: legacy appended, then
% the nested + SS vent build without the CL resonance, then with it. If clvoct
% moves level_c materially, the resonance is not the free improvement the
% contrast/maperr trade suggested (+52% chi for 0.5% maperr).

cfg = { 'm=3b stock (CONTROL)  ', 3, 0
        'm=4 appended (legacy) ', 4, 1
        'm=4 nested+SS, no res ', 4, 2
        'm=4 default (clvoct 4)', 4, 3 };

fprintf('\n  target level_c 5 (=1.62 %%/dB), band [3.5 6.5] (=1.26-1.88 %%/dB)\n');
fprintf('  SS=0.10 measured 36.75 (=3.67 %%/dB) -- the failure that motivated this\n');
fprintf('  CHECK ROW 1 FIRST: the fitted scaffolding should be IN BAND.\n\n');
fprintf('  config                 | level_c | %%/dB  | slope | maperr | in band\n');
fprintf('%s\n', repmat('-',1,68));
R = struct('nm',{},'lc',{},'pdb',{},'maperr',{});
for i = 1:size(cfg,1)
    if (cfg{i,2} == 3)
        pa = modpar26(3); pa.chsz = [0.95 0.05 1.00];
    else
        pa = modpar26(4);
        switch cfg{i,3}
            case 1      % legacy appended
                pa.chsz = [0.95 0.05 1.00 0.05]; pa.nested = 0;
                for fn = {'clvent','clvoct','clvtgt'}
                    if (isfield(pa,fn{1})), pa = rmfield(pa,fn{1}); end
                end
            case 2      % nested + SS vent, resonance disabled
                pa = rmfield(pa,'clvoct');
        end
    end
    lc = NaN; sl = NaN; pdb = NaN; mp = NaN;
    try
        evalc('m = abr_metric(pa, false);');
        lc = m.level_c; sl = m.slope;
        if (isfinite(lc) && lc > 0), pdb = 100*(lc^(1/100)-1); end
    catch e
        fprintf('  %-22s | abr_metric FAILED: %s\n', cfg{i,1}, e.message);
    end
    try
        S = score26(pa, 'fast', false); mp = S.maperr;
    catch
    end
    inb = 'no';
    if (isfinite(lc) && lc>=3.5 && lc<=6.5), inb = 'YES'; end
    fprintf('  %-22s | %7.2f | %5.2f | %5.3f | %6.1f | %s\n', ...
        cfg{i,1}, lc, pdb, sl, mp, inb);
    R(end+1).nm=cfg{i,1}; R(end).lc=lc; R(end).pdb=pdb; R(end).maperr=mp; %#ok<SAGROW>
    save('m4_levelc.mat','R');
end
fprintf(['\n  If row 1 is in band and rows 3-4 are not, m=4''s level dependence is\n' ...
         '  broken and the amplifier result stands on an incomplete metric set --\n' ...
         '  the same way the SS=0.10 result did. If rows 3-4 are IN band, m=4\n' ...
         '  survives a check that SS=0.10 failed, which is the first evidence the\n' ...
         '  nested construction is sound on an axis it was never tuned for.\n']);
disp('M4_LEVELC_DONE');
