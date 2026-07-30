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

% ---- 2. third DOF made INERT must reproduce the 2-DOF model -------------
% REDEFINED 2026-07-29 for the lab frame. The old form froze the third DOF with
% m5o*1e6 and that is no longer the inert limit: with dref(3)=false, a3 = s3/m5,
% so m5 -> inf pins d3 IN THE LAB (d3 -> 0) and k5 becomes a spring to GROUND,
% which stiffens the BM instead of decoupling. Measured: the old form returns
% 1.00038, a 100% difference, and it is not a regression -- the test was wrong.
%
% The lab-frame inert limit is to REMOVE THE ATTACHMENT: k5 -> 0 and r5 -> 0.
% Then spas = 0, its reaction on d1 vanishes, d3 is driven by nothing and stays
% at zero, so the model must reproduce plain 2-DOF m=3 exactly.
% Expected value is still ~0, so the stored baseline stays valid.
v = NaN;
try
    p3  = modpar26(3); p3.d3int = 0;              % legacy 2-DOF
    p3b = modpar26(3);                            % 3-DOF, attachment removed
    p3b.k5o = 0; p3b.r5o = 0;
    v = wdiff(p3, p3b);
catch e, R.err2 = e.message; end
[R,npass,nfail] = gchk(R,B,'g2_m3_m3b',v,0,1e-9,'k5=r5=0 -> 2DOF','~0 identity',verbose,npass,nfail);

% ---- 3. third DOF LIVE must actually change the model -------------------
% REPLACED 2026-07-29. The old g3 tested m=3b -> m=4 at clcouple=0, an
% equivalence that NO LONGER EXISTS BY DESIGN: the m<4 correction (lab frame +
% reaction) was kept while m=4 was reverted, so the two force laws now differ on
% purpose. Its construction was also unsound -- clcouple=0 zeroes d3's pressure
% pickup while chamber 4 still exists, leaving an UNCOUPLED fluid compartment
% that diverges at sample 3. Repairing it would mean asserting something false.
%
% Replaced by the INVERSE of g2, and by the check whose absence cost this whole
% project. For its entire history the third DOF at m<4 was EXACTLY INERT: d3 live
% and unfrozen gave rel-to-2DOF = 0.000e+00, bit-identical, so every m=3b result
% ever recorded was a 2-DOF result and the old g2 "passed" because there was
% nothing to freeze. Nothing watched for that.
%
% g2 and g3 are now a proper pair: INERT when it should be (k5=r5=0 -> identity),
% LIVE when it should be (default k5/r5 -> a real difference). One alone is not
% enough; the old suite had only the first kind.
% NEW KEY NAME on purpose: comparing to the old g3_m3b_m4 baseline (3.5784e-10)
% would be a false comparison, since the quantity means the opposite now.
v = NaN;
try
    p3  = modpar26(3); p3.d3int = 0;              % 2-DOF reference
    p3b = modpar26(3);                            % 3-DOF, attachment LIVE
    v = wdiff(p3, p3b);
catch e, R.err3 = e.message; end
R.g3_dof3_live = v;
LIVEFLOOR = 1e-3;                                  % below this, d3 does nothing
ok3 = isfinite(v) && v > LIVEFLOOR;
if (verbose)
    st = 'FAIL'; if (ok3), st = 'PASS'; end
    fprintf('  %-29s | %15.6g | %-13s | %s (live)\n', ...
            'DOF3 live rel-to-2DOF', v, sprintf('>%.0e', LIVEFLOOR), st);
end
if (ok3), npass = npass + 1; else, nfail = nfail + 1; end

% ---- 4. scored champions -----------------------------------------------
% exp4 is the FALLBACK used only when no baseline file exists; gchk prefers the
% baseline whenever it is present. m3b was 499.3 until 2026-07-29, when fdm26
% gained the third DOF at m=3 -- before that fdm26 scored a 2-DOF model while
% tdm26 ran a 3-DOF one, so 499.3 measured a model nobody was running. Updated so
% a capture-from-scratch does not report a spurious failure, and so the printed
% "expected" column stops disagreeing with the value actually compared against.
exp4 = [104.6 1263.95 525.2];
lbl4 = {'g4_maperr_m1','g4_maperr_m3b','g4_maperr_m4'};
lbl8 = {'g8_amp_m1','g8_amp_m3b','g8_amp_m4'};
ampv = NaN(1,3);
for k = 1:3
    v = NaN;
    try
        if (k==1), pk = modpar26(1); elseif (k==2), pk = modpar26(3); else, pk = modpar26(4); end
        S = score26(pk, 'fast', false); v = S.maperr;
        % AMPLIFIER GAIN, captured from the SAME score26 call -- it was already
        % being computed and thrown away, so this costs nothing.
        if (isfield(S,'amp_gain') && isnumeric(S.amp_gain) && ~isempty(S.amp_gain))
            ampv(k) = double(S.amp_gain(1));
        end
    catch e, R.(sprintf('err4_%d',k)) = e.message; end
    [R,npass,nfail] = gchk(R,B,lbl4{k},v,exp4(k),0.15,lbl4{k},sprintf('%.1f',exp4(k)),verbose,npass,nfail);
end

% ---- 8. AMPLIFIER MUST EXIST -------------------------------------------
% ADDED 2026-07-29. This gate exists because of a specific failure: the m=4
% "nested, no vent" configuration scores maperr 353.5 (the BEST unfitted number
% in the project) with a healthy 4.14-octave map and maxRe 0.0 (perfectly
% sub-critical) -- and amp_gain = -0.01 dB. It has NO COCHLEAR AMPLIFIER. Active
% and passive responses are identical.
%
% Every other column rewards that. maperr improves 33%, stability looks perfect,
% and tip-tail contrast goes UP (16.07 vs 6.41). Deleting the mechanism the model
% exists to study is invisible to, and rewarded by, every gate above. The same
% signature is on record for the clvoct=0.5 prescription ("leaves no amplifier at
% all, amp +0.05"). A near-zero amp_gain with excellent stability is this
% model's characteristic failure, and it reads as success.
%
% So: an ABSOLUTE FLOOR, not just a baseline comparison. A baseline captured on a
% broken configuration would enshrine amp_gain ~ 0 as correct; the floor cannot.
% Measured references: m=4 default +56.59, m=3b fitted +81.15.
AMPFLOOR = 1.0;                        % dB; below this there is no amplifier
for k = 1:3
    a = ampv(k);
    R.(lbl8{k}) = a;
    ok = isfinite(a) && a > AMPFLOOR;
    ref = NaN; if (isfield(B,lbl8{k}) && isfinite(B.(lbl8{k}))), ref = B.(lbl8{k}); end
    % also fail a >20% drift from a baseline that itself cleared the floor
    if (ok && isfinite(ref) && ref > AMPFLOOR && abs(a-ref) > 0.2*abs(ref)), ok = false; end
    if (verbose)
        st = 'FAIL'; if (ok), st = 'PASS'; end
        note = sprintf('>%.1f dB', AMPFLOOR);
        if (isfinite(ref)), note = sprintf('%s, ~%.2f', note, ref); end
        fprintf('  %-29s | %15.4f | %-13s | %s (amp)\n', ...
                sprintf('amplifier %s', strrep(lbl8{k},'g8_amp_','')), a, note, st);
    end
    if (ok), npass = npass + 1; else, nfail = nfail + 1; end
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

% ---- 7. STABILITY: every configuration must stay finite ------------------
% ADDED 2026-07-29 after m=4 diverged in tdm26 for a whole session while this
% gate reported PASS. score26 'fast' runs FDM26, which is frequency-domain and
% has NO stability limit, so a maperr can look healthy while the time march
% blows up. Nothing here watched the march itself. Now something does.
% Reported as max|d1| per config: NaN/Inf or a huge value both fail.
STAB = {'m=1',1,{}; 'm=2',2,{}; 'm=3',3,{'d3int',0}; 'm=3b',3,{}; 'm=4',4,{}};
ISV7 = [1136 1005 840 655 466 273 80];
for k = 1:size(STAB,1)
    v = NaN;
    try
        pk = modpar26(STAB{k,2}); ex = STAB{k,3};
        for j = 1:2:numel(ex), pk.(ex{j}) = ex{j+1}; end
        pk.isv = ISV7;
        evalc('Sk = tdm26(0,pk,0,0);'); dk = Sk.d1(:);
        if (all(isfinite(dk))), v = max(abs(dk)); else, v = Inf; end
    catch e, R.(sprintf('err7_%d',k)) = e.message; end
    % Pass band is generous: this catches DIVERGENCE, not drift. Anything above
    % 1e-1 m of BM displacement is nonphysical by orders of magnitude.
    ok = isfinite(v) && v < 1e-1;
    R.(sprintf('g7_finite_%s', strrep(STAB{k,1},'=',''))) = v;
    if (verbose)
        st = 'FAIL'; if (ok), st = 'PASS'; end
        fprintf('  %-29s | %15.6g | finite <1e-1  | %s (stab)\n', ...
                sprintf('stability %s max|d1|',STAB{k,1}), v, st);
    end
    if (ok), npass = npass + 1; else, nfail = nfail + 1; end
end

% ---- summary / capture --------------------------------------------------
% SURFACE CAPTURED EXCEPTIONS. Every gate above stores its error text in an
% err* field and NOTHING ever printed them, so a thrown exception and a
% divergent waveform both surfaced as a bare NaN. Distinguishing those two cost
% most of 2026-07-29: I read NaN as divergence, then "corrected" myself to
% exception, then back, with the answer sitting unprinted in R the whole time.
if (verbose)
    ef = fieldnames(R); ef = ef(strncmp(ef,'err',3));
    if (~isempty(ef))
        fprintf('\n  CAPTURED EXCEPTIONS (a NaN gate above is explained here):\n');
        for k = 1:numel(ef)
            msg = R.(ef{k}); msg = strrep(msg, sprintf('\n'), ' ');
            fprintf('    %-8s %s\n', ef{k}, msg(1:min(100,end)));
        end
    else
        fprintf('\n  no exceptions captured (any NaN above is DIVERGENCE, not a throw)\n');
    end
end

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
% SAY WHY, do not just return NaN. A silent NaN here is indistinguishable from a
% thrown exception in the caller's try/catch, which is exactly how m=4's
% divergence hid behind "g3 failed" for a full session on 2026-07-29.
if (numel(a)~=numel(b))
    fprintf('    [wdiff] LENGTH MISMATCH %d vs %d\n', numel(a), numel(b)); d = NaN; return;
end
if (~all(isfinite(a)) || ~all(isfinite(b)))
    ia = find(~isfinite(a),1); ib = find(~isfinite(b),1);
    fprintf('    [wdiff] DIVERGED, not an exception: A first at %s, B first at %s\n', ...
            num2str(ia), num2str(ib)); d = NaN; return;
end
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
