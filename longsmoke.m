% LONGSMOKE -- 3-evaluation dry run of longfit's exact configuration.
%
% Cheap insurance on an 18-hour commitment. The one thing that must be true
% before launching: parfit26 accepts the STRUCT warm start and actually starts
% at ace=+1.00, rather than silently falling back to modpar26(1). That fallback
% is real and has bitten this project already -- sweep1b passed no 'warm', hit
% the chamber-count mismatch at parfit26.m:157, and started from the baseline
% without saying so. If it fires here the run would spend 18 h re-deriving a
% result we already have.
%
% PASS = the start line reads slope ~0.4415 and maperr ~107.7 (the ace=+1.00
% point), NOT slope 0.663 / maperr 104.6 (the modpar26(1) baseline).
RUNDIR='/private/tmp';
L=load('sweep_nch1b.mat'); pa=L.R.pa; pa.ace=1.00;
fprintf('\n== longfit smoke: expect start slope ~0.4415, maperr ~107.7 ==\n');
fprintf('   (a fallback to modpar26(1) would read slope 0.663, maperr 104.6)\n\n');
R = parfit26(1, struct('maxfe',3, 'wgain',0.01, 'warm',pa, 'verbterm',true, ...
                       'out',fullfile(RUNDIR,'longsmoke.mat')));
fprintf('\n   ace in returned pa: %+.4f   (want ~+1.00, NOT -0.40)\n', R.pa.ace);
fprintf('   ace at pv index 20: %+.4f\n', R.pv(20));
disp('LONGSMOKE_DONE');
