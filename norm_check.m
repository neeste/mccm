% Did removing the fixed-sum chsz normalization (a) leave m=1/2/3 untouched, and
% (b) make the m=4 reduction gate exact?
ISV=[1136 1005 840 655 466 273 80];
fprintf('\n(a) REGRESSION -- m=1/2/3 must be unchanged (their chsz already sums to 2)\n');
for m=[1 2 3]
    pa=modpar26(m); pa.isv=ISV;
    pl=pa; pl.chsznorm=1;                       % legacy normalization
    evalc('Sa=tdm26(0,pa,0,0);'); evalc('Sb=tdm26(0,pl,0,0);');
    fprintf('   m=%d  max|d1 new - d1 legacy| = %.3e  (scale %.3e)\n', ...
            m, max(abs(Sa.d1(:)-Sb.d1(:))), max(abs(Sb.d1(:))));
end
fprintf('\n(b) GATE -- m=4 at clcouple=0 vs m=3b\n');
p3b=modpar26(3); p3b.isv=ISV;
p40=modpar26(4); p40.clcouple=0; p40.isv=ISV;
evalc('S3b=tdm26(0,p3b,0,0);'); evalc('S40=tdm26(0,p40,0,0);');
d1d=max(abs(S3b.d1(:)-S40.d1(:))); d2d=max(abs(S3b.d2(:)-S40.d2(:)));
d3d=max(abs(S3b.d3(:)-S40.d3(:))); pdd=max(abs(S3b.ped(:)-S40.ped(:)));
fprintf('   max|d1 diff| = %.3e   (was 4.389e-04, scale %.3e)\n', d1d, max(abs(S3b.d1(:))));
fprintf('   max|d2 diff| = %.3e   (was 6.649e-04)\n', d2d);
fprintf('   max|d3 diff| = %.3e   (was 6.572e-04)\n', d3d);
fprintf('   max|ped diff|= %.3e\n', pdd);
if (max([d1d d2d d3d pdd]) < 1e-18)
    fprintf('\n   GATE PASSED: m=4 now reduces to m=3b EXACTLY.\n');
else
    fprintf('\n   still differs -- normalization was not the only cause.\n');
end
disp('NORM_CHECK_DONE');
