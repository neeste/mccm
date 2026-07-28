p2=modpar26(2); p3=modpar26(3); s2=p3; s2.m=2; s3=p2; s3.m=3; s3.chsz=[1 1 1];
cfg={{'native2',p2,0.445},{'native3',p3,0.77},{'swap2',s2,0.530},{'swap3',s3,0.578}};
fprintf('\n%-9s %7s %9s %8s %6s %6s   %s\n','config','d','expected','hiCFcon','nval','R2','BF sequence');
for c=1:numel(cfg)
    m=tiptail_metric(cfg{c}{2});
    if (isempty(m.BF)), fprintf('%-9s FAILED: %s\n',cfg{c}{1},m.msg); continue; end
    fprintf('%-9s %7.3f %9.3f %8.1f %6d %6.2f   %s\n', ...
            cfg{c}{1}, m.d, cfg{c}{3}, m.chi, m.nvalid, m.r2, num2str(m.BF,'%7.2f'));
    if (~isempty(m.msg)), fprintf('           note: %s\n', m.msg); end
end
disp('TTV2_DONE');
