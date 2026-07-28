m=checkcode('parfit26.m','-struct');
syn=m(arrayfun(@(z)any(cellfun(@(s)~isempty(strfind(lower(z.message),s)),{'parse','syntax','unbalanced'})),m));
if (~isempty(syn))
    for k=1:numel(syn), fprintf('SYNTAX L%d %s\n',syn(k).line,syn(k).message); end
    return
end
fprintf('checkcode OK\n\n');
% Exercise the NEW surface objective on the current best parameters (1 eval).
o.warm='parfit26_recip.mat'; o.maxfe=1; o.surface=1; o.wsurf=1;
o.wshoulder=0; o.out='tmp_surface_check.mat';
R=parfit26(3,o);
fprintf('\n--- surface-objective check ---\n');
fprintf('Delta      = %.2f ms  -> predicted neural delay a = %.2f ms (2013 assumed 5.00)\n', R.delta, 5+R.delta);
fprintf('band-RMS   = %.3f ms over %d valid cells\n', R.surf_rms, R.surf_nvalid);
fprintf('slope/level (diagnostic only): d=%.3f  level_c=%.2f\n', R.mf.slope, R.mf.level_c);
fprintf('python point-fit predicted Delta in [0.58 (2013), 0.83 (1988)]; band loss should sit in/near that range\n');
disp('VERIFY_SURFACE_DONE');
