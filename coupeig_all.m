% Print maxRe (ALL modes) vs maxRe_osc for the m3form=1 cases that diverged.
% If maxRe > 0 while maxRe_osc < 0, the instability is STATIC and maxRe_osc misses it.
fprintf('\n  m3form gam   maxRe(all)   maxRe_osc    top-mode f(kHz)   time-march\n');
for gam=[0 0.30 0.70]
  for mf=[0 1]
    pa=modpar26(3); pa.gam=gam; pa.m3form=mf;
    evalc('E=tdm26(''coupeig'',struct(''pa'',pa));');
    ftop=abs(imag(E.lam(1)))/2/pi/1000;
    % quick divergence check
    pa2=pa; pa2.isv=1391:-10:11;
    dv='finite'; try, evalc('S=tdm26(0,pa2,0,0);'); if(any(~isfinite(S.d1(:)))),dv='DIVERGES';end; catch, dv='THREW'; end
    fprintf('  %d     %.2f   %+10.1f  %+10.1f   %10.2f      %s\n', mf,gam,E.maxRe,E.maxRe_osc,ftop,dv);
  end
end
disp('COUPEIG_ALL_DONE');
