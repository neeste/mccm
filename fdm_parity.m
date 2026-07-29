function R = fdm_parity(verbose)
%FDM_PARITY  Do fdm26 and tdm26 describe the SAME partition topology?
%
% The capstone design note requires "fdm/tdm parity checked at the partition
% boundary". This is that check, and it tests the right invariant.
%
% WHAT IS COMPARED: the TOPOLOGY -- which chamber pair each DOF spans, and with
% which orientation. That is what must never diverge. When fdm26 lacked the
% nested/clvent port, maperr described the appended model while tdm26 amplified
% the nested one, and every m=4 map number was meaningless for a whole session.
% This test exists so that cannot recur silently.
%
% WHAT IS NOT COMPARED: the numerical values. The two solvers pose different
% equations -- fdm26 builds Y = T*(zk\D) in an impedance formulation, tdm26
% builds s = s_int - Df*p with the mass ratios applied separately -- so the
% scalings legitimately differ:
%     m=2      fdm [1 -1]      tdm [0.5 -0.5]     factor of two
%     m=4 d3   fdm [0 -1 0 1]  tdm [0 -cc 0 cc]   clcouple is tdm-side
%     B        fdm [-1;0;1]    tdm [1;0;-1]       sign convention
% Demanding numerical equality would make this test permanently red for reasons
% that are not faults. Demanding nothing would make it vacuous.
%
% Row 3 of fdm26's m>=3 D is its volume-conservation slot (all zeros), which has
% no tdm26 counterpart; it is excluded rather than counted as a mismatch.

if (nargin<1), verbose = true; end
cfg = { 'm=1',                        modpar26(1)
        'm=2',                        modpar26(2)
        'm=3',                        modpar26(3)
        'm=4 appended',               setf(modpar26(4),'nested',0)
        'm=4 nested',                 setf(modpar26(4),'nested',1)
        'm=4 nested, vent->SS',       setf(setf(modpar26(4),'nested',1),'clvtgt',2) };
R = struct('name',{},'ok',{},'note',{});
np = 0; nf = 0;
if (verbose)
    fprintf('\n  FDM/TDM PARITY at the partition boundary\n');
    fprintf('  Comparing TOPOLOGY (which chambers each DOF spans), not scaling.\n\n');
    fprintf('  configuration            | dofs | topology | note\n');
    fprintf('  %s\n', repmat('-',1,66));
end
for i = 1:size(cfg,1)
    nm = cfg{i,1}; pa = cfg{i,2}; note = '';
    try
        F = fdm26(struct('macroD',1,'pa',pa));
        cp = fake_cp(pa);
        C  = macro_couple(pa, cp);
        nd = min(size(F.topo,1), C.ndof);
        a = sign(F.topo(1:nd,:));
        b = sign(C.Df(1:nd,1:size(F.topo,2)));
        keep = any(a~=0,2) | any(b~=0,2);      % drop all-zero rows (conservation slot)
        ok = isequal(a(keep,:), b(keep,:));
        if (~ok)
            note = 'TOPOLOGY DIVERGES';
            if (verbose)
                fprintf('  %-24s | %4d | MISMATCH | see below\n', nm, nd);
                disp('      fdm topo:'); disp(a(keep,:));
                disp('      tdm topo:'); disp(b(keep,:));
            end
        end
    catch e
        ok = false; nd = -1; note = ['THREW: ' e.message];
    end
    if (ok), np=np+1; st='match'; else, nf=nf+1; st='DIFFER'; end
    if (verbose && ok), fprintf('  %-24s | %4d | %-8s | %s\n', nm, nd, st, note); end
    R(end+1).name=nm; R(end).ok=ok; R(end).note=note; %#ok<AGROW>
end
if (verbose)
    fprintf('\n  %d match, %d differ\n', np, nf);
    if (nf==0)
        fprintf('  PARITY HOLDS -- both solvers describe the same partition topology.\n');
    else
        fprintf('  PARITY BROKEN -- the solvers describe DIFFERENT models. Any maperr\n');
        fprintf('  compared against a tdm26 amplifier result is meaningless until fixed.\n');
    end
end
end

function p = setf(p, f, v), p.(f) = v; end

function cp = fake_cp(pa)
n = pa.n; x = linspace(0,1,n).';
cp.m1 = 0.0075*exp(-0.20*x); cp.m2 = 0.0360*exp(-0.08*x);
if (pa.m>=4 || (isfield(pa,'d3int') && pa.d3int)), cp.m5 = 0.0360*exp(-0.08*x); end
if (isfield(pa,'clvent') && pa.clvent>0 && isfield(cp,'m5'))
    cp.clvm = cp.m5/pa.clvent;
    kv = 0; if (isfield(pa,'clvk')), kv = pa.clvk; end
    cp.clvk = kv*ones(n,1); cp.clvr = zeros(n,1);
end
end
