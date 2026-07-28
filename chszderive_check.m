% Does chszderive actually HOLD sum(chsz)=2 through a search?
%
% Testing the constraint, not that the option is accepted. The wstat smoke test
% passed earlier for a reason unrelated to wstat (maxRe never crossed statol),
% so here the fit is given the chsz indices to move and the sum is checked after.
%
% ARM A  chsz free (chszderive off) -- the sum is EXPECTED to drift off 2.00.
%        If it does not drift, the fit is not moving chsz and the test is inert.
% ARM B  chszderive=3 (SV derived)  -- the sum MUST be 2.00 to machine precision
%        no matter where the search went.
%
% m=4, pure map path, small budget. fitidx includes 31,32,34 (ST,SS,CL); index
% 33 (SV) is included deliberately in arm A and should be auto-dropped in arm B.
pa = modpar26(4); pa.clvent = 0.5;
b = struct('wslope',0,'wlevel',0,'wlcb',0,'skipabr',1,'cheapstab',1, ...
           'maptol',0,'fitidx',[31 32 33 34],'warm',pa,'maxfe',25, ...
           'pin',struct('clvent',0.5),'out','chszderive_check.mat');
for k = 1:2
    o = b;
    if (k==2), o.chszderive = 3; end
    lbl = {'A chsz FREE','B SV DERIVED'};
    fprintf('\n  ===== %s =====\n', lbl{k});
    R = parfit26(4, o);
    cz = R.pa.chsz;
    fprintf('  chsz  = [%s]\n', strtrim(sprintf('%.4f ', cz)));
    fprintf('  sum   = %.10f   (deviation %.2e)\n', sum(cz), abs(sum(cz)-2));
    fprintf('  maperr= %.1f\n', R.Rf.maperr);
end
fprintf(['\n  ARM A should show a sum that has DRIFTED off 2.00 (else the fit is\n' ...
         '  not moving chsz and this proves nothing). ARM B must be 2.00 to\n' ...
         '  machine precision regardless of where the search went.\n']);
disp('CHSZDERIVE_CHECK_DONE');
