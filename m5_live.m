for mm=[1 8]
  pa=modpar26(4); pa.chsz(4)=0; pa.m5o=0.0360276*mm;
  fprintf('\nset pa.m5o = %.6g  (x%g)\n', pa.m5o, mm);
  evalc('S=tdm26(0,pa,0,0);');
  if (isfield(S,'cp'))
    c=S.cp;
    if (isfield(c,'m5'))
      mu3=c.m1./max(c.m5,1e-12);
      fprintf('  cp.m5  min=%.6g max=%.6g\n', min(c.m5), max(c.m5));
      fprintf('  mu3    min=%.6g max=%.6g\n', min(mu3), max(mu3));
    else, fprintf('  cp HAS NO FIELD m5\n'); end
    fprintf('  S.pa.m5o = %.6g   (what the model actually used)\n', S.pa.m5o);
    fprintf('  S.pa.chsz = [%s]\n', num2str(S.pa.chsz,'%6.3f'));
  else, fprintf('  S has no cp field\n'); end
end
disp('M5_LIVE_DONE');
