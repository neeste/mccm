% REDUCTION GATE FOR THE REBUILT m=4.
%
% LADDER
%   m=3   3 chambers, 2 DOFs
%   m=3b  3 chambers, 3 DOFs, d3 INTERNAL          (gate PASSED: exact identity)
%   m=4   4 chambers, 3 DOFs, d3 gets CL           (this gate)
%
% pa.clcouple in [0,1] scales d3's FLUID coupling. At 0 the CL chamber is
% decoupled and d3 is internal, so m=4 must reproduce m=3b EXACTLY. That is the
% gate. chsz(4)->0 does NOT serve: there mu3 cancels algebraically and d3 still
% carries a pressure constraint.
%
% The BM force sign is now FIXED BY THE GATE rather than chosen: m=3b has -act
% on the BM, so the rebuilt m=4 carries the same pair. The legacy m=4 had it
% inverted, which the bisection showed was why its amplifier was dead.
%
% Beyond the gate, sweeping clcouple up shows what giving d3 a compartment
% actually does, one step at a time, which no previous m=4 experiment could see.

ISV = [1136 1005 840 655 466 273 80];
p3b = modpar26(3);  p3b.isv = ISV;              % m=3b (d3int on by default now)
evalc('S3b = tdm26(0,p3b,0,0);');

fprintf('\n  === GATE: m=4 at clcouple=0 must equal m=3b ===\n');
p40 = modpar26(4); p40.clcouple = 0; p40.isv = ISV;
evalc('S40 = tdm26(0,p40,0,0);');
d1d = max(abs(S3b.d1(:)-S40.d1(:)));
d2d = max(abs(S3b.d2(:)-S40.d2(:)));
d3d = max(abs(S3b.d3(:)-S40.d3(:)));
pdd = max(abs(S3b.ped(:)-S40.ped(:)));
fprintf('   max|d1 diff| = %.3e   (scale %.3e)\n', d1d, max(abs(S3b.d1(:))));
fprintf('   max|d2 diff| = %.3e\n', d2d);
fprintf('   max|d3 diff| = %.3e\n', d3d);
fprintf('   max|ped diff|= %.3e\n', pdd);
if (max([d1d d2d d3d pdd]) < 1e-18)
    fprintf('   GATE PASSED: the rebuilt m=4 reduces to m=3b exactly.\n');
else
    fprintf('   *** GATE FAILED -- the reduction is not exact. ***\n');
end

fprintf('\n  === what giving d3 a compartment actually does ===\n');
fprintf('  clcouple | amp d1  amp d2  amp d3 | maxRe     | range mono | maperr\n');
fprintf('%s\n', repmat('-',1,74));
for cc = [0 0.05 0.2 0.5 1.0]
    pa = modpar26(4); pa.clcouple = cc;
    try
        S = score26(pa, 'fast', false);
        fprintf('  %8.2f | %+6.2f %+7.2f %+7.2f | %+9.1f | %5.2f %-4s | %6.1f\n', ...
            cc, S.amp_gain, S.amp_d2, S.amp_d3, S.maxRe, S.bf_range, S.bf_mono, S.maperr);
    catch e
        fprintf('  %8.2f | FAILED: %s\n', cc, e.message);
    end
end
fprintf(['\n  m=3 reference: amp d1 +81.15  d2 +84.17 | maxRe +19.3 | maperr 499.3\n' ...
         '  Watch where amp and stability depart from m=3b as the compartment is\n' ...
         '  introduced. That is the step the bisection said adds uncalibrated gain.\n']);
disp('M4_GATE_DONE');
