% THREE-COORDINATE AMPLIFIER GAIN. Regression + the m=4 question.
% CONTROL: m=2 native d1 gain must still read +39.11 dB (the validated value).
% If it does, the d1 path is untouched and the new d2/d3 columns can be read.
% m=1/m=2 have no separate micromechanical DOFs so their gain must appear in d1;
% for m=4 the question is whether gain sits in d2/d3 while d1 stays small.
cfg={ {'m=2 native (control)', modpar26(2)}, ...
      {'m=3 native',           modpar26(3)} };
b2=modpar26(4).m2o; br1=modpar26(4).r1o;
p4=modpar26(4); p4.m2o=b2*32;                       cfg{end+1}={'m=4 m2o x32', p4};
p4b=p4; p4b.r1o=br1*0.5;                            cfg{end+1}={'m=4 +r1 x0.5', p4b};
p4c=p4b; p4c.ohcsgn=-1; p4c.ohcgain=0.01;           cfg{end+1}={'m=4 sgn-1 og.01', p4c};
fprintf('\n  config             | d1(BM)  | d2(shear)| d3(OC ht)| max     | maxRe\n');
fprintf('%s\n',repmat('-',1,74));
for c=1:numel(cfg)
    nm=cfg{c}{1}; pa=cfg{c}{2};
    try
        S=score26(pa,'fast',false);
        fprintf('  %-18s | %+7.2f | %+8.2f | %+8.2f | %+7.2f | %+8.1f\n', ...
            nm, S.amp_gain, S.amp_d2, S.amp_d3, S.amp_max, S.maxRe);
    catch e
        fprintf('  %-18s | FAILED: %s\n', nm, e.message);
    end
end
fprintf('\n  CONTROL: m=2 d1 must be +39.11 dB. If it moved, the d1 path broke.\n');
fprintf('  m=4 READ: large d2/d3 with small d1 => the amplifier WORKS locally and\n');
fprintf('  the problem is BM coupling. Small everywhere => the amplifier is broken.\n');
disp('AMP3_TEST_DONE');
