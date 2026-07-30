function R = macro_shadow(verbose)
% MACRO_SHADOW  Stage 1 gate: does macro_couple reproduce the hand-written stamps?
%
% Builds Df/Dq/mu/B via macro_couple and, INDEPENDENTLY, evaluates the literal
% formulas transcribed from tdm26's xpnd_q (721-766) and fold_p (768-815), then
% asserts they agree on random inputs. tdm26 itself is UNTOUCHED -- this proves
% the matrix form before Stage 2 switches anything over.
%
% Random s and p rather than a real march: a real run exercises one trajectory,
% random inputs exercise the whole linear map. If the two agree on random data
% they agree everywhere, and the test cannot pass by hitting a lucky state.
%
% Coverage is every branch tdm26 can take:
%   m=1, m=2, m=3, m=3b (d3int)
%   m=4 appended, clcouple = 1 and 0
%   m=4 nested, vent to SS and SV, with and without the resonant clvk
%   m=4 nested + clcouple ~= 1  -- the known non-reciprocal case, checked to
%                                 FAIL reciprocity but still MATCH the code
%
% A note on that last one: the transcription below reproduces tdm26's asymmetry
% deliberately (fold applies cc under nested, xpnd does not). If a future edit
% makes tdm26 reciprocal, THIS TEST WILL FAIL -- which is correct. It is pinning
% current behaviour, not endorsing it.

if (nargin<1), verbose = true; end
rng(7);                                  % deterministic: a flaky gate is useless
cfg = shadow_cfgs();
R = struct('name',{},'dfold',{},'dxpnd',{},'recip',{});
np = 0; nf = 0;

if (verbose)
    fprintf('\n  MACRO SHADOW -- does macro_couple reproduce the hand-written stamps?\n');
    fprintf('  tdm26 is NOT modified. Random s,p over every branch.\n\n');
    fprintf('  configuration                  | fold err  | xpnd err  | recip | status\n');
    fprintf('  %s\n', repmat('-',1,72));
end

for i = 1:numel(cfg)
    c = cfg(i);
    try
        pa = c.pa; cp = fake_cp(pa);
        C  = macro_couple(pa, cp);
        n  = pa.n; nch = C.nch; ndof = C.ndof;
        s  = randn(n, ndof);             % internal partition forces
        p  = randn(n*nch, 1);            % chamber pressures, interleaved by place
        qst= randn;

        [fa, qa] = via_matrix(C, s, p, qst, n);
        [fb, qb] = via_literal(pa, cp, s, p, qst, n, nch, ndof);

        df = relerr(fa, fb); dq = relerr(qa, qb);
        ok = (df < 1e-12) && (dq < 1e-12);
    catch e
        df = NaN; dq = NaN; ok = false; C.reciprocal = false;
        if (verbose), fprintf('  %-30s | THREW: %s\n', c.name, e.message); end
    end
    if (ok), st = 'PASS'; np = np+1; else, st = 'FAIL'; nf = nf+1; end
    if (verbose)
        rc = 'no '; if (C.reciprocal), rc = 'yes'; end
        fprintf('  %-30s | %9.2e | %9.2e | %-5s | %s\n', c.name, df, dq, rc, st);
    end
    R(end+1).name = c.name; R(end).dfold = df; R(end).dxpnd = dq; R(end).recip = C.reciprocal; %#ok<AGROW>
end

% ---- FLUID BLOCK: is a2 a SECOND copy of the same topology? --------------
% ADDED 2026-07-29. The test above checks the partition<->fluid TRANSFER
% (fold_p / xpnd_q). It never checked the fluid operator's own chamber-to-chamber
% block, which macro26 writes out BY HAND. That block is redundant with
% macro_couple:
%
%     a2(k) == diag(L_p + L_c) + Dq*diag(mu(k,:))*Df
%
% Two sources of truth for one topology is why pa.rlsplit is ill-posed: it
% changes Df/Dq and the hand-written block does not follow, over-determining the
% system (diverges at sample 26, invariant under a 40x change in CL area).
%
% This uses the REAL macro26 cp, not fake_cp, because only macro26 builds a2.
[npf, nff] = deal(0);
if (verbose)
    fprintf('\n  FLUID BLOCK -- does a2 equal diag(L) + Dq*diag(mu)*Df?\n');
    fprintf('  config          | max abs err | verdict\n');
    fprintf('  %s\n', repmat('-',1,52));
end
for mm = [3 4]
    lbl = sprintf('modpar26(%d)', mm); e = NaN;
    try
        pa = modpar26(mm); cp = macro26(pa); C = cp.mc;
        nn = pa.n; nc = C.nch; dx = pa.xl/(nn-1);
        af = 1; if (isfield(pa,'aflom_fac')), af = pa.aflom_fac; end
        afl = cp.ac / (af * pa.rho * dx);
        e = 0;
        for k = round(0.3*nn):round(0.1*nn):round(0.7*nn)   % interior places only
            L = zeros(1,nc);
            for c = 1:nc, L(c) = (afl(k-1)+afl(k))*pa.chsz(c)./cp.abmom(k); end
            A = reshape(cp.a2(k,:), nc, nc) - diag(L);
            P = C.Dq * diag(C.mu(k,:)) * C.Df;
            e = max(e, max(abs(A(:)-P(:))));
        end
    catch ex
        e = NaN;
    end
    ok = isfinite(e) && e < 1e-12;
    if (ok), npf = npf+1; else, nff = nff+1; end
    if (verbose)
        v = 'DIFFERS'; if (ok), v = 'identical'; end
        fprintf('  %-15s | %11.3e | %s\n', lbl, e, v);
    end
end
if (verbose)
    fprintf('\n  identical => the hand-written block can be REPLACED by the product,\n');
    fprintf('  which is what would make pa.rlsplit well-posed. DIFFERS at m=4 would\n');
    fprintf('  localise the known nested+clcouple non-reciprocity.\n');
end

if (verbose)
    fprintf('\n  %d passed, %d FAILED  (transfer)\n', np, nf);
    fprintf('  %d passed, %d FAILED  (fluid block)\n', npf, nff);
    if (nf==0)
        fprintf('  Stage 1 GATE HOLDS -- the matrix form reproduces every branch.\n');
        fprintf('  Safe to proceed to Stage 2 (switch xpnd_q/fold_p over).\n');
    else
        fprintf('  Stage 1 GATE BROKEN -- macro_couple does not match tdm26. Do not proceed.\n');
    end
end
end

% =========================================================================
function [fs, q] = via_matrix(C, s, p, qst, n)
% The proposed form: s - Df*p, and Dq*(mu.*s) + qst*B.
P = reshape(p, C.nch, n).';                       % (n x nch), place-major
fs = s - P * C.Df.';                              % (n x ndof)
ws = C.mu .* s;                                   % mass-weighted forces
Q  = ws * C.Dq.';                                 % (n x nch)
Q(1,:) = Q(1,:) + qst * C.B.';                    % stapes enters at the base only
q  = reshape(Q.', n*C.nch, 1);
end

function [fs, q] = via_literal(pa, cp, s, p, qst, n, nch, ndof)
% Literal transcription of tdm26's xpnd_q and fold_p. Deliberately verbose and
% branchy -- it is the thing being reproduced, so it must not be tidied.
m = pa.m;
nest = isfield(pa,'nested') && pa.nested;
cc = 1; if (isfield(pa,'clcouple')), cc = pa.clcouple; end
s1 = s(:,1); s2 = s(:,2);
s3 = zeros(n,1); if (ndof>=3), s3 = s(:,3); end
s4 = zeros(n,1); if (ndof>=4), s4 = s(:,4); end
mu2 = cp.m1 ./ max(cp.m2,1e-12);
mu3 = zeros(n,1); if (isfield(cp,'m5')), mu3 = cp.m1 ./ max(cp.m5,1e-12); end
muv = zeros(n,1); if (isfield(cp,'clvm')), muv = cp.m1 ./ max(cp.clvm,1e-12); end
qx  = zeros(n*nch,1); qs = [qst; zeros(n-1,1)];
fs  = s;

if (m<=1)
    fs(:,1) = s1 - p;
    qx = s1 + qs;
elseif (m==2)
    j2=(1:n)*2; j1=j2-1;
    fs(:,1) = s1 - (p(j1) - p(j2))/2;
    qx(j1) = s1 + qs;
elseif (m==3)
    j3=(1:n)*3; j2=j3-1; j1=j2-1;
    fs(:,1) = s1 - (p(j1) - p(j3));
    fs(:,2) = s2 - (p(j3) - p(j2));
    qx(j1) =  s1 + qs;
    qx(j2) = -s2 .* mu2;
    qx(j3) = -s1 + s2 .* mu2 - qs;
else
    j4=(1:n)*4; j3=j4-1; j2=j4-2; j1=j4-3;
    if (nest), fs(:,1) = s1 - (p(j1) - p(j4));
    else,      fs(:,1) = s1 - (p(j1) - p(j3)); end
    fs(:,2) = s2 - (p(j3) - p(j2));
    fs(:,3) = s3 - cc*(p(j4) - p(j2));           % cc in BOTH branches
    if (nest)
        qx(j1) =  s1 + qs;
        qx(j2) = -s2 .* mu2 - s3 .* mu3;         % bare mu3 -- no cc (the asymmetry)
        qx(j3) =        s2 .* mu2 - qs;
        qx(j4) = -s1 + s3 .* mu3;                % bare mu3 -- no cc
        if (ndof>=4)
            vt = 3; if (isfield(pa,'clvtgt')), vt = pa.clvtgt; end
            jt = j1; if (vt==2), jt=j2; elseif (vt==3), jt=j3; end
            fs(:,4) = s4 - (p(j4) - p(jt));
            qx(j4) = qx(j4) + s4 .* muv;
            qx(jt) = qx(jt) - s4 .* muv;
        end
    else
        m3q = cc .* mu3;
        qx(j1) =  s1 + qs;
        qx(j2) = -s2 .* mu2 - s3 .* m3q;
        qx(j3) = -s1 + s2 .* mu2 - qs;
        qx(j4) =  s3 .* m3q;
    end
end
q = qx;
end

function cp = fake_cp(pa)
% Minimal cp with the fields the interface reads. Smooth place-dependent masses
% so mu varies along the cochlea, as it does in the real model -- a constant mu
% would hide an indexing error.
n = pa.n; x = linspace(0,1,n).';
cp.m1 = 0.0075 * exp(-0.20*x);
cp.m2 = 0.0360 * exp(-0.08*x);
if (pa.m>=4 || (isfield(pa,'d3int') && pa.d3int))
    cp.m5 = 0.0360 * exp(-0.08*x);
end
if (isfield(pa,'clvent') && pa.clvent>0 && isfield(cp,'m5'))
    cp.clvm = cp.m5 / pa.clvent;
    kv = 0; if (isfield(pa,'clvk')), kv = pa.clvk; end
    cp.clvk = kv * ones(n,1);
    cp.clvr = zeros(n,1);
end
end

function d = relerr(a, b)
sc = max(max(abs(a(:))), max(abs(b(:)))); if (sc<=0), sc = 1; end
d = max(abs(a(:)-b(:))) / sc;
end

function cfg = shadow_cfgs()
% Every branch tdm26 can take. n is kept small: this is a linear-algebra
% identity check, not a physics run, so 41 places exercises the same code.
N = 41;
cfg = struct('name',{},'pa',{});
mk = @(nm,pa) struct('name',nm,'pa',pa);

p = base_pa(1,N);                      cfg(end+1) = mk('m=1', p);
p = base_pa(2,N);                      cfg(end+1) = mk('m=2 (symmetry constraint)', p);
p = base_pa(3,N); p.d3int = 0;         cfg(end+1) = mk('m=3 (2 DOF)', p);
p = base_pa(3,N); p.d3int = 1;         cfg(end+1) = mk('m=3b (internal d3)', p);
p = base_pa(4,N);                      cfg(end+1) = mk('m=4 appended, cc=1', p);
p = base_pa(4,N); p.clcouple = 0;      cfg(end+1) = mk('m=4 appended, cc=0 (gate)', p);
p = base_pa(4,N); p.clcouple = 0.35;   cfg(end+1) = mk('m=4 appended, cc=0.35', p);
p = base_pa(4,N); p.nested = 1; p.clvtgt = 2; p.clvent = 3;
                                       cfg(end+1) = mk('m=4 nested, vent->SS', p);
p = base_pa(4,N); p.nested = 1; p.clvtgt = 3; p.clvent = 0.5;
                                       cfg(end+1) = mk('m=4 nested, vent->SV', p);
p = base_pa(4,N); p.nested = 1; p.clvtgt = 2; p.clvent = 3; p.clvk = 1e7;
                                       cfg(end+1) = mk('m=4 nested, resonant vent', p);
p = base_pa(4,N); p.nested = 1; p.clvtgt = 2; p.clvent = 3; p.clcouple = 0.5;
                                       cfg(end+1) = mk('m=4 nested + cc=0.5 (non-recip)', p);
end

function pa = base_pa(m, n)
pa = struct('m',m,'n',n);
if (m>=4), pa.d3int = 0; end
end
