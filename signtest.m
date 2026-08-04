% SIGNTEST -- does the OHC RC pole INJECT or ABSORB energy at the BM?
% Direct net-work measurement, which is what the question actually is. The
% earlier impedance argument (rcsign) was invalid: it read Re(-gam*za) as if the
% active term were a self-impedance, but it is a CROSS-coupling (row-1
% coefficient of V2), and this model already amplifies +41 dB with Re exactly 0.
%
% ANCHOR FIRST: with the pole OFF this model amplifies +41 dB, so its ohcBM must
% be POSITIVE. If it is not, the sign derivation in micro26 is wrong and every
% number below is meaningless.
%
% PREDICTION UNDER TEST: the pole is supposed to add a second 90 deg ABOVE its
% corner, turning negative stiffness into negative damping. If so, ohcBM should
% RISE for f > fc and be roughly unchanged for f << fc. If it FALLS above the
% corner, the pole is adding dissipation and the mechanism runs backwards.
PROJ = '/Users/neely/Library/CloudStorage/OneDrive-FatherFlanagan''sBoysHome/STN/mod/mccm';
addpath(PROJ); cd(PROJ);
L=load('fit_nch3_surface.mat'); pa=L.R.pa;
if (isfield(pa,'hbmode')), pa=rmfield(pa,'hbmode'); end
fc=1000; tau=1/(2*pi*fc);
fprintf('\n  corner fc = %d Hz.  ohcBM > 0 = energy INTO the BM.\n', fc);
fprintf('\n  %8s | %14s %14s | %10s\n','f (Hz)','ohcBM off','ohcBM on','ratio on/off');
for f=[250 500 1000 2000 4000 8000]
    r=struct('fr',f/1000,'lv',60,'pa',pa);
    S0=tdm26('wnr1',r,0,0);
    p=pa; p.ohctau=tau; r1=struct('fr',f/1000,'lv',60,'pa',p);
    S1=tdm26('wnr1',r1,0,0);
    a=S0.dgn.ohcBM; b=S1.dgn.ohcBM;
    fprintf('  %8.0f | %14.4e %14.4e | %10.3f\n', f, a, b, b/a);
end
fprintf('\n  ANCHOR: the OFF column must be POSITIVE at every frequency (this\n');
fprintf('  model amplifies +41 dB). If it is not, micro26''s sign is wrong.\n');
fprintf('  VERDICT: ratio > 1 above 1000 Hz = the pole INJECTS (mechanism works).\n');
fprintf('           ratio < 1 above 1000 Hz = the pole ABSORBS (runs backwards).\n');
disp('SIGNTEST_DONE');
