% IS m=3b BIT-IDENTICAL TO m=3 IN d1/d2?
% d3 in m=3b has NO back-action: s1 and s2 depend only on d1,d2,v1,v2; d3 is
% internal so it enters no fluid equation; a(i1) and a(i2) never read it. So d3
% is a SHADOW coordinate -- it shows what the OC height would do under m=3
% dynamics without perturbing them. If that is right the reduction gate is an
% EXACT IDENTITY, which is the strongest possible test of the new plumbing:
% any difference in d1/d2 is a bug in the m=3b implementation, full stop.
ISV=[1136 1005 840 655 466 273 80];
p3=modpar26(3);   p3.isv=ISV;
pb=modpar26c3b;   pb.isv=ISV;
evalc('S3=tdm26(0,p3,0,0);');
evalc('Sb=tdm26(0,pb,0,0);');
d1d=max(abs(S3.d1(:)-Sb.d1(:)));
d2d=max(abs(S3.d2(:)-Sb.d2(:)));
pd =max(abs(S3.ped(:)-Sb.ped(:)));
fprintf('\n  max |d1(m=3) - d1(m=3b)|  = %.3e\n', d1d);
fprintf('  max |d2(m=3) - d2(m=3b)|  = %.3e\n', d2d);
fprintf('  max |ped(m=3) - ped(m=3b)|= %.3e\n', pd);
fprintf('  scale: max|d1| = %.3e\n', max(abs(S3.d1(:))));
if (d1d==0 && d2d==0 && pd==0)
    fprintf('\n  IDENTICAL. d3 is a pure shadow coordinate; the plumbing is correct.\n');
elseif (d1d < 1e-18)
    fprintf('\n  identical to roundoff. plumbing correct.\n');
else
    fprintf('\n  *** DIFFERS -- d3 is feeding back somewhere it should not. BUG. ***\n');
end
if (isfield(Sb,'d3') && any(Sb.d3(:)~=0))
    fprintf('  m=3b d3 (shadow RL) max = %.3e  vs d1 %.3e  (ratio %.2f)\n', ...
            max(abs(Sb.d3(:))), max(abs(Sb.d1(:))), max(abs(Sb.d3(:)))/max(abs(Sb.d1(:))));
else
    fprintf('  NOTE: d3 is all zeros -- it is not being driven at all.\n');
end
disp('M3B_IDENT_DONE');
