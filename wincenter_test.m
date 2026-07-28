function wincenter_test
% Where does the 2 kHz coherent reflection form?  Fix the stimulus at 2 kHz
% (WNR~2.68, 2xWNR~5.4 ms) and move a +-0.75-oct roughness band across CF place.
% The center that maximizes emag AND gives OAE~5.4 ms is the reflection region.
fr=2;
fprintf('\n=== window-center sweep, stimulus %g kHz (2xWNR~5.4 ms), rough_amp=0.1 ===\n',fr);
fprintf(' win_fc(kHz)  WNR    OAE    emag(dB)\n');
for fc=[0.5 1 2 3 4 6]
  pr.fr=fr; pr.lv=60; pr.pa=modpar26(3); pr.pa.oae=1; pr.pa.rough_amp=1e-1;
  pr.pa.rough_fc=fc; pr.pa.rough_foct=0.75;
  S=tdm26('wnr1',pr,0,0);
  fprintf('  %6.2f    %5.2f  %5.2f   %6.1f\n', fc, S.tpk, S.oae, S.oam);
end
disp('WINCENTER_DONE');
