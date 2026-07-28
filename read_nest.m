L=load('nesting_test.mat'); R=L.R;
fprintf('\n  config                     | maperr | range mono | contr | amp    | maxRe\n');
fprintf('%s\n',repmat('-',1,82));
for i=1:numel(R)
    S=R(i).S;
    fprintf('  %-26s | %8.1f | %5.2f %-4s | %5.1f | %+6.2f | %+7.1f\n', ...
        R(i).tag, S.maperr, S.bf_range, S.bf_mono, S.contrast, S.amp_gain, S.maxRe);
end
fprintf('\n  (amp NaN where the passive reference gam=0 diverges -- known for m3form=1)\n');
disp('READ_NEST_DONE');
