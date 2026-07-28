S2=score26(modpar26(2),'fast',false); p4=modpar26(4); p4.m2o=p4.m2o*32;
S4=score26(p4,'fast',false);
fprintf('\nFOLD CALIBRATION (magnitude, octaves):\n');
fprintf('  m=2 native (gold standard): range %.2f  fold %.3f -> %s\n', S2.bf_range,S2.bf_fold,S2.bf_mono);
fprintf('  m=4 m2o x32               : range %.2f  fold %.3f -> %s\n', S4.bf_range,S4.bf_fold,S4.bf_mono);
disp('SMOKE2_DONE');
