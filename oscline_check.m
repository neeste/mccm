% Does the corrected stability line catch the case that misled me?
% refit_m4_map.mat holds the VOID result: maxRe +27.6 (unstable) with
% maxRe_osc -2.3. The OLD line printed only the osc value and labelled it
% "sub-critical", so reading parfit26's own summary would have passed that run.
% Warm-starting from those parameters with maxfe=1 reproduces the exact case.
% EXPECT: "maxRe : +27.6 (UNSTABLE) [osc -2.3] <== osc alone would say
% sub-critical". Anything reporting sub-critical means the fix did not take.
L = load('refit_m4_map.mat');
pa = L.R.pa;
fprintf('\n  loaded void-run pa: n=%d chsz=[%s] clvent=%.2f\n', ...
    pa.n, strtrim(sprintf('%.2f ', pa.chsz)), pa.clvent);
o = struct('wslope',0,'wlevel',0,'wlcb',0,'skipabr',1,'cheapstab',0, ...
           'maptol',0,'fitidx',[1 2],'warm',pa,'maxfe',1, ...
           'pin',struct('chsz',pa.chsz,'clvent',pa.clvent), ...
           'out','oscline_check.mat');
parfit26(4, o);
fprintf(['\n  The stability line above must show the FULL maxRe and call it\n' ...
         '  UNSTABLE. The old line showed only osc (-2.3) and said sub-critical.\n']);
disp('OSCLINE_CHECK_DONE');
