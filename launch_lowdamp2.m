% Continuation of the low-damping fit to bring maperr under the 185 tolerance.
% Why it stalled at 191.3: with wmap=0.01 the 6.3-point excess contributed only
% 0.063 to J=0.994 (~1.5%), so the optimizer had no incentive. wmap is a HINGE
% (wmap*max(0,maperr-185)), so once maperr<185 the term and its gradient vanish
% -- a large weight forces the crossing without distorting the fit beyond it.
L=load('parfit26_lowdamp.mat'); pa=L.R.pa;      % continue from the current best
fprintf('continuing from: maperr=%.1f  bandRMS=%.3f\n', L.R.Rf.maperr, L.R.surf_rms);
o.warm     = pa;
o.hbmode   = 'bm';
o.surface  = 1;  o.wsurf = 1;
o.maptol   = 185; o.wmap = 0.30;          % 30x stronger; hinge => no effect once under
o.wcrit    = 0.05;
o.wshoulder= 0;
o.maxfe    = 120;
o.out      = 'parfit26_lowdamp2.mat';
R = parfit26(3,o);
fprintf('\n=== CONTINUATION: maperr=%.1f (t<=185)  bandRMS=%.3f ms  Delta=%.2f  osc=%+.1f ===\n', ...
        R.Rf.maperr, R.surf_rms, R.delta, R.S.maxRe_osc);
fprintf('    maperr under tolerance: %d   (was 191.3)\n', R.Rf.maperr<=185);
fprintf('    bandRMS vs previous 0.931: %+.3f\n', R.surf_rms-0.931);
f=R.mf.f(:); lat=R.mf.lat; i2=find(abs(f-2)<0.01,1);
fprintf('2 kHz latency steps: ');
for b=1:size(lat,2)-1, fprintf('%+.1f%% ',100*(lat(i2,b+1)-lat(i2,b))/lat(i2,b)); end
fprintf('\n');
disp('LOWDAMP2_DONE');
