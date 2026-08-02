% GPGUARD -- prove the gampro grade cannot be reached by accident.
%
% The specific hazard, caught during review rather than by a test: at nch=1/2 the
% grade lives at pv index 33, and the DEFAULT fitidx already carries 33 as a
% chsz(3) that those models do not have. If the out-of-range clip were widened to
% admit the grade, every default nch=1/2 call would have started fitting gampro
% silently, and the 90.97 result would no longer be reproducible from its own
% recorded settings. So: default must still DROP 33 and fit 15 params; only
% opts.fitgampro may add it, taking the count to 16.
fprintf('\n== gampro reachability guard ==\n');
tmp = tempname;
for nch = [1 3]
    for fg = [false true]
        s = evalc(sprintf('R=parfit26(%d,struct(''maxfe'',1,''fitgampro'',%d,''out'',''%s''));', ...
                          nch, fg, [tmp sprintf('_%d_%d.mat',nch,fg)]));
        nfit  = regexp(s, 'fitting (\d+) params', 'tokens', 'once');
        nfit  = str2double(nfit{1});
        saidgp = ~isempty(strfind(s, 'fitting gampro grade')); %#ok<STREMP>
        sawdrop= ~isempty(strfind(s, 'dropping'));             %#ok<STREMP>
        fprintf('  nch=%d fitgampro=%d | params %2d | gampro announced %d | clip fired %d\n', ...
                nch, fg, nfit, saidgp, sawdrop);
    end
end
fprintf('\n  EXPECT nch=1: 15 params / gampro 0 / clip 1   then  16 / 1 / 1\n');
fprintf('  EXPECT nch=3: 16 params / gampro 0 / clip 0   then  17 / 1 / 0\n');
disp('GPGUARD_DONE');
