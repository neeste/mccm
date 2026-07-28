% Re-fit from REDUCED DAMPING under the surface objective, with maperr enforced.
% Rationale: r1x0.5 gave the best latency-step distribution seen (-25.7/-24.0/
% -54.9 vs baseline -26.7/-6.8/-50.9) but wrecked tuning (maperr 1991.9). r1o is
% itself a fitted parameter, so scaling it alone breaks maperr because nothing
% compensates. This asks whether the improved distribution survives once the
% other parameters can re-optimize under the maperr constraint.
L=load('parfit26_recip.mat'); pa=L.R.pa;
pa.r1o = pa.r1o*0.5;                      % the lever that moved the distribution
o.warm     = pa;
o.hbmode   = 'bm';                        % BM-peak metric (never the default drive)
o.surface  = 1;  o.wsurf = 1;             % 16-cell band objective, free Delta
o.maptol   = 185; o.wmap = 0.01;          % strong: maperr must come down from ~1992
o.wcrit    = 0.05;                        % keep sub-critical
o.wshoulder= 0;
o.maxfe    = 150;
o.out      = 'parfit26_lowdamp.mat';
R = parfit26(3,o);
fprintf('\n=== LOW-DAMP RE-FIT: maperr=%.1f (t<=185)  bandRMS=%.3f ms  Delta=%.2f  osc=%+.1f ===\n', ...
        R.Rf.maperr, R.surf_rms, R.delta, R.S.maxRe_osc);
f=R.mf.f(:); lat=R.mf.lat; i2=find(abs(f-2)<0.01,1);
fprintf('2 kHz latency steps: ');
for b=1:size(lat,2)-1, fprintf('%+.1f%% ',100*(lat(i2,b+1)-lat(i2,b))/lat(i2,b)); end
fprintf(' (target -28%% each)\n');
disp('LOWDAMP_DONE');
