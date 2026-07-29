function m3dof_check
%M3DOF_CHECK  The m<3 third DOF: what changed, and what must NOT have.
%
% SN: "when m<3 the micro-mechanics should mimic m=4 via interaction between SS
% and CL without longitudinal coupling."
%
% This is the FIRST behaviour-changing step of the refactor, so the check has two
% halves that must not be conflated:
%   (A) REGRESSION -- every existing configuration must be untouched. m=1/m=2 at
%       the default dof=2 still take the legacy law, so arch_gate must read 9/9.
%       If it does not, the change leaked into configurations it had no business
%       touching and the result is void regardless of (B).
%   (B) NEW CAPABILITY -- m=1 with dof=3 must run, and must differ from m=1 with
%       dof=2. A new configuration that runs but returns the OLD answer would
%       mean the third DOF is inert, which is the failure most likely to look
%       like success.

fprintf('\n  (A) REGRESSION -- existing configurations must be untouched\n');
R = arch_gate(false);
fprintf('      arch_gate: %d passed, %d FAILED\n', R.npass, R.nfail);
if (R.nfail > 0)
    fprintf('      *** LEAKED into existing configurations -- (B) is moot ***\n');
    fns = fieldnames(R);
    for i=1:numel(fns)
        if (strncmp(fns{i},'g',1)), fprintf('        %-16s %.6g\n', fns{i}, R.(fns{i})); end
    end
    return
end
fprintf('      m=1/m=2 at dof=2 still take the legacy law: confirmed\n');

fprintf('\n  (B) NEW CAPABILITY -- m=1 with dof=3\n');
p2 = modpar26(1); p2.isv = [1136 1005 840 655 466 273 80];
p3 = p2; p3.dof = 3;
fprintf('      k5o present in modpar26(1) now: %d\n', isfield(p2,'k5o'));
try
    evalc('S2 = tdm26(0,p2,0,0);');
    evalc('S3 = tdm26(0,p3,0,0);');
    a2 = max(abs(S2.d1(:))); a3 = max(abs(S3.d1(:)));
    fin = all(isfinite(S3.d1(:)));
    has3 = isfield(S3,'d3') && any(S3.d3(:)~=0);
    rel = max(abs(S3.d1(:)-S2.d1(:))) / max(a2,eps);
    fprintf('      dof=2  max|d1| %.4e\n', a2);
    fprintf('      dof=3  max|d1| %.4e   finite %d   d3 active %d\n', a3, fin, has3);
    fprintf('      relative d1 difference: %.3e\n', rel);
    if (~fin)
        fprintf('      RUNS BUT DIVERGES -- not usable\n');
    elseif (rel < 1e-12)
        fprintf('      *** THIRD DOF IS INERT -- identical to dof=2, no effect ***\n');
    elseif (~has3)
        fprintf('      *** d3 state never moves -- the DOF exists but carries nothing ***\n');
    else
        fprintf('      WORKS: runs, finite, d3 active, and the answer CHANGED.\n');
    end
catch e
    fprintf('      FAILED: %s\n', e.message);
end
fprintf(['\n  Note m=1 dof=3 is NOT expected to match anything previously recorded --\n' ...
         '  it is a configuration that has never existed. The point of (A) is that\n' ...
         '  making it reachable did not disturb what already worked.\n']);
end
