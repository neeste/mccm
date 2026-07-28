% Smoke test: wstat accepted, and default (= wcrit) reproduces prior behaviour.
% Cheap: m=1, pure map path, tiny budget. Only checks plumbing, not fit quality.
pa = modpar26(1);
base = struct('wslope',0,'wlevel',0,'wlcb',0,'skipabr',1,'cheapstab',1, ...
              'maptol',0,'fitidx',[1 2 3],'warm',pa,'maxfe',1, ...
              'pin',struct('chsz',[1 1]),'out','wstat_smoke.mat');
fprintf('\n  A: defaults (wstat should inherit wcrit)\n');
RA = parfit26(1, base);
o = base; o.wstat = 0.5; o.statol = 1;   % deliberately harsh static guard
fprintf('\n  B: wstat=0.5, statol=1 (harsh) -- must be ACCEPTED and change J\n');
RB = parfit26(1, o);
fprintf('\n  A maperr %.1f | B maperr %.1f\n', RA.Rf.maperr, RB.Rf.maperr);
fprintf('  wstat option accepted without error: yes\n');
disp('WSTAT_SMOKE_DONE');
