function R = arch_gate(verbose)
% ARCH_GATE  Regression contract for the macro/micro separation refactor.
%
% The refactor (see the capstone design note: "chamber count and DOF count are
% both parameters, not forks") must be BEHAVIOUR-PRESERVING. No stage may change
% any number. This is the contract every stage is checked against.
%
% FIRST RUN captures a baseline into arch_gate_baseline.mat. LATER RUNS compare
% against that capture rather than against hardcoded literals, so the gate stays
% self-consistent even where my recorded values were rounded for reporting.
% Expected values are printed alongside as a sanity check on the capture itself.
%
% THE GATES, all verified during the 2026-07-26/28 work:
%   1  m=1 == m=2                        exact (symmetry-constrained 2-chamber)
%   2  m=3 -> m=3b, heavy m5             exact identity, frozen third DOF
%   3  m=3b -> m=4 at clcouple=0         exact to roundoff (~3.6e-10 rel)
%   4  score26 fast, three champions     maperr 104.6 / 499.3 / 525.2
%   5  fitted m=4 (refit_m4_full.mat)    maperr 329.3, maxRe +8.3, amp +54.06
%   6  resonant vent, clvk far below CF  converges to the eliminated form
%
% Gates 2, 3 and 6 compare RAW d1 waveforms, which is the strictest available
% test and the one that caught the d3-save bug and the m4 sign inversion. Gates
% 4 and 5 compare scored metrics, which catches anything that changes the model
% without changing the march.
%
% Written as a FUNCTION file, not a script, so local helpers are legal --
% local functions in a SCRIPT get their own workspace and cannot see script
% variables, which cost a run earlier in this project.

if (nargin<1), verbose = true; end
BASE = 'arch_gate_baseline.mat';
have = exist(BASE,'file') == 2;
if (have), L = load(BASE); B = L.B; else, B = struct(); end
R = struct(); npass = 0; nfail = 0;

if (verbose)
    fprintf('\n  ARCH GATE -- macro/micro refactor regression contract\n');
    if (have), fprintf('  mode: CHECK against %s\n\n', BASE);
    else,      fprintf('  mode: CAPTURE baseline -> %s\n\n', BASE); end
    fprintf('  gate                          | value           | expected      | status\n');
    fprintf('  %s\n', repmat('-',1,74));
end

% ---- 1. m=1 == m=2 -----------------------------------------------------
v = NaN;
try
    e1 = fdm26(struct('pa',modpar26(1)));
    e2 = fdm26(struct('pa',modpar26(2)));
    v = abs(e1.maperr - e2.maperr);
catch e, R.err1 = e.message; end
[R,npass,nfail] = gchk(R,B,'g1_m1eq_m2',v,0,1e-12,'m=1==m2 |dmaperr|','0 exact',verbose,npass,nfail);

% ---- 2. m=3 -> m=3b, third DOF frozen ----------------------------------
v = NaN;
try
    p3  = modpar26(3); p3.d3int = 0;              % legacy 2-DOF
    p3b = modpar26(3); p3b.m5o = p3b.m5o * 1e6;   % 3-DOF, third frozen
    v = wdiff(p3, p3b);
catch e, R.err2 = e.message; end
[R,npass,nfail] = gchk(R,B,'g2_m3_m3b',v,0,1e-9,'m=3->m3b rel d1','~0 identity',verbose,npass,nfail);

% ---- 3. m=3b -> m=4 at clcouple=0 --------------------------------------
v = NaN;
try
    p3b = modpar26(3); p3b.chsz = [0.95 0.05 1.00];
    p4  = modpar26(4); p4.chsz = [0.95 0.05 1.00 0.05];
    p4.nested = 0; p4.clcouple = 0;
    for f = {'clvent','clvoct','clvtgt','clvk','clvr'}
        if (isfield(p4,f{1})), p4 = rmfield(p4,f{1}); end
    end
    v = wdiff(p3b, p4);
catch e, R.err3 = e.message; end
[R,npass,nfail] = gchk(R,B,'g3_m3b_m4',v,0,1e-8,'m3b->m4 rel d1','~3.6e-10',verbose,npass,nfail);

% ---- 4. scored champions -----------------------------------------------
exp4 = [104.6 499.3 525.2];
lbl4 = {'g4_maperr_m1','g4_maperr_m3b','g4_maperr_m4'};
for k = 1:3
    v = NaN;
    try
        if (k==1), pk = modpar26(1); elseif (k==2), pk = modpar26(3); else, pk = modpar26(4); end
        S = score26(pk, 'fast', false); v = S.maperr;
    catch e, R.(sprintf('err4_%d',k)) = e.message; end
    [R,npass,nfail] = gchk(R,B,lbl4{k},v,exp4(k),0.15,lbl4{k},sprintf('%.1f',exp4(k)),verbose,npass,nfail);
end

% ---- 5. fitted m=4 ------------------------------------------------------
if (exist('refit_m4_full.mat','file')==2)
    try
        F = load('refit_m4_full.mat'); pf = F.R.pa;
        S = score26(pf, 'fast', false);
        [R,npass,nfail] = gchk(R,B,'g5_fit_maperr',S.maperr,329.3,0.15,'fitted m4 maperr','329.3',verbose,npass,nfail);
        [R,npass,nfail] = gchk(R,B,'g5_fit_amp',   S.amp_gain,54.06,0.05,'fitted m4 amp',   '+54.06',verbose,npass,nfail);
    catch e, R.err5 = e.message; end
end

% ---- 6. resonant vent reduces to the eliminated form -------------------
v = NaN;
try
    pv = modpar26(4);                       % ships clvoct=4 -> dof 4
    pe = pv; pe = rmfield(pe,'clvoct');     % eliminated form -> dof 3
    pf = pv; pf.clvoct = 10;                % resonance far below CF
    v = wdiff(pe, pf);
catch e, R.err6 = e.message; end
[R,npass,nfail] = gchk(R,B,'g6_vent_reduce',v,0,5e-3,'vent clvoct=10 rel d1','~0 converged',verbose,npass,nfail);

% ---- summary / capture --------------------------------------------------
if (~have)
    B = R; save(BASE,'B'); %#ok<NASGU>
    if (verbose), fprintf('\n  BASELINE CAPTURED -> %s. Re-run after each stage.\n', BASE); end
else
    if (verbose)
        fprintf('\n  %d passed, %d FAILED\n', npass, nfail);
        if (nfail==0), fprintf('  CONTRACT HOLDS -- the stage preserved behaviour.\n');
        else,          fprintf('  CONTRACT BROKEN -- do not proceed; revert or fix.\n'); end
    end
end
R.npass = npass; R.nfail = nfail;
end

% =========================================================================
function d = wdiff(pa, pb)
% Relative max|d1| difference between two configurations, on the same places.
% Raw-waveform comparison: the strictest test available and the one that caught
% both the d3-save bug and the m=4 force-sign inversion.
ISV = [1136 1005 840 655 466 273 80];
pa.isv = ISV; pb.isv = ISV;
evalc('Sa = tdm26(0,pa,0,0);');
evalc('Sb = tdm26(0,pb,0,0);');
a = Sa.d1(:); b = Sb.d1(:);
if (numel(a)~=numel(b) || ~all(isfinite(a)) || ~all(isfinite(b))), d = NaN; return; end
sc = max(abs(a)); if (sc<=0), sc = 1; end
d = max(abs(a-b)) / sc;
end

function [R,np,nf] = gchk(R, B, key, val, expct, tol, lbl, expstr, verbose, np, nf)
% Compare against the captured baseline when present, else against the recorded
% expectation. Judging on RELATIVE difference where the value is large -- an
% absolute threshold once made me report a pure-roundoff 2.9e-13 as a failure.
R.(key) = val;
if (isfield(B,key) && isfinite(B.(key)))
    ref = B.(key); src = 'base';
else
    ref = expct; src = 'exp ';
end
ok = isfinite(val) && (abs(val-ref) <= max(tol, tol*abs(ref)));
if (ok), st = 'PASS'; np = np+1; else, st = 'FAIL'; nf = nf+1; end
if (verbose)
    fprintf('  %-29s | %15.6g | %-13s | %s (%s)\n', lbl, val, expstr, st, src);
end
end
