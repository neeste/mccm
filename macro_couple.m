function C = macro_couple(pa, cp)
% MACRO_COUPLE  The macromechanics/micromechanics interface, made explicit.
%
% The macro/micro separation (capstone design note: "chamber count and DOF count
% are both parameters, not forks"). Built once per run by macro26 and carried as
% cp.mc; tdm26's xpnd_q and fold_p consume it. macro_shadow.m proved it
% reproduces every hand-written branch it replaced (11 configurations, errors
% 0 to 6.8e-17).
%
% THE INTERFACE. tdm26's xpnd_q and fold_p are hand-written expansions of two
% matrix operations:
%
%   fold_p :  s = s_internal - Df*p                 pressure differences drive DOFs
%   xpnd_q :  q(:,c) = sum_k Dq(c,k)*mu(:,k).*s(:,k) + qst*B(c)
%
% Df is (ndof x nch), Dq is (nch x ndof), B is (nch x 1). mu is (n x ndof) --
% PLACE-DEPENDENT, since the mass ratios m1/m2 and m1/m5 vary along the cochlea,
% which is why mu is carried separately rather than folded into Dq.
%
% fdm26.m:916-931 already carries this as D and B; this brings tdm26 to the same
% structure rather than inventing one.
%
% TWO PLACES WHERE Dq IS *NOT* Df' -- both faithfully preserved here:
%
% (1) m=2 IS DELIBERATELY ASYMMETRIC. fold takes HALF the pressure difference
%     (Df = [0.5 -0.5]) while xpnd injects into chamber 1 only (Dq = [1;0]), and
%     tdm26.m:547 replaces the second chamber equation with a constraint row.
%     That is the symmetry constraint SN described: "the 1-chamber model actually
%     represents 2-chamber physics with the added constraint of symmetry". Not a
%     bug; do not "fix" it. It is what makes m=1 == m=2 hold exactly.
%
% (2) NESTED + clcouple IS NON-RECIPROCAL, and this one IS a latent bug.
%     fold_p applies ccf to the d3 pickup in BOTH branches (tdm26.m:798-799),
%     but xpnd_q applies ccq only in the non-nested branch -- its nested branch
%     uses bare mu3 (tdm26.m:744,746). So with nested=1 and clcouple~=1 the
%     pickup is scaled and the injection is not, and the operator loses
%     reciprocity. It has never bitten because every nested configuration runs
%     clcouple=1 (the default), and clcouple is only used as the m=3b->m=4
%     reduction knob, which is a non-nested path.
%     PRESERVED HERE so Stage 1 stays behaviour-preserving. C.reciprocal reports
%     whether Df' and Dq agree, so a caller can assert it. Fixing it is a
%     separate, visible decision -- see the note at the end of this file.
%
% Returns C with fields:
%   nch, ndof   counts
%   Df          (ndof x nch)  pressure-difference selector
%   Dq          (nch x ndof)  volume-injection topology (mu applied separately)
%   mu          (n x ndof)    mass ratios m1/m_k, place-dependent
%   B           (nch x 1)     stapes injection pattern
%   dref        (ndof x 1)    logical: does DOF k's acceleration add a(i1)?
%   reciprocal  logical: does Dq == Df' hold for this configuration?

m = pa.m; n = pa.n;
nest = isfield(pa,'nested') && pa.nested;
cc   = 1; if (isfield(pa,'clcouple')), cc = pa.clcouple; end
has3 = (m>=4) || (isfield(pa,'d3int') && pa.d3int) || (isfield(pa,'dof') && pa.dof>=3);
has4 = (m>=4) && isfield(cp,'clvm') && ~isempty(cp.clvm) && ...
       isfield(cp,'clvk') && any(cp.clvk(:)~=0);

% ---- DOF count and reference convention --------------------------------
% dref: DOF k's acceleration is a(ik) = a(i1) + s_k/M_k when true, else s_k/M_k.
% This is a MICROMECHANICS property that used to be keyed on CHAMBER COUNT --
% exactly the entanglement the refactor exists to break. It now follows the
% micromechanics VARIANT instead; see the note below.
ndof = 2 + double(has3) + double(has4);
dref = false(ndof,1);
% The reference convention follows the MICROMECHANICS VARIANT, not the chamber
% count: whenever the m=3/m=4 force law runs, d2 and d3 are ABSOLUTE (measured
% from the BM). m<3 with dof=2 keeps the legacy RELATIVE d2 untouched, which is
% what preserves m=1 == m=2.
use3 = (m>=3) || (isfield(pa,'dof') && pa.dof>=3);
if (use3), dref(2) = true; if (has3), dref(3) = true; end, end
% dref(4) stays false: the vent state is a fluid FLOW, independent of BM motion.

% ---- mass ratios, place-dependent --------------------------------------
mu = ones(n, ndof);                                   % mu(:,1) == 1 by definition
if (ndof>=2 && isfield(cp,'m2')), mu(:,2) = cp.m1 ./ max(cp.m2, 1e-12); end
if (ndof>=3 && isfield(cp,'m5')), mu(:,3) = cp.m1 ./ max(cp.m5, 1e-12); end
if (ndof>=4 && isfield(cp,'clvm')), mu(:,4) = cp.m1 ./ max(cp.clvm, 1e-12); end

% ---- topology per chamber count ----------------------------------------
nch = max(m,1);
Df = zeros(ndof, nch); Dq = zeros(nch, ndof); B = zeros(nch,1);
switch true
    case (m<=1)                                  % 1 chamber, d2 internal
        Df(1,1) = 1;  Dq(1,1) = 1;  B(1) = 1;
    case (m==2)                                  % see note (1): asymmetric
        Df(1,:) = [0.5 -0.5];
        Dq(1,1) = 1;                             % chamber 2 carries a constraint
        B(1) = 1;
    case (m==3)
        Df(1,:) = [ 1  0 -1];                    % d1: P_ST - P_SV
        Df(2,:) = [ 0 -1  1];                    % d2: P_SV - P_SS
        % d3 (m=3b) is INTERNAL: no fluid compartment, so Df(3,:) stays zero.
        Dq(1,1) =  1;
        Dq(2,2) = -1;
        Dq(3,1) = -1;  Dq(3,2) = 1;
        B = [1;0;-1];
    otherwise                                    % m>=4
        if (nest), Df(1,:) = [1 0 0 -1];         % d1: P_ST - P_CL
        else,      Df(1,:) = [1 0 -1 0];         % d1: P_ST - P_SV
        end
        Df(2,:) = [0 -1  1  0];                  % d2: P_SV - P_SS
        Df(3,:) = [0 -cc 0 cc];                  % d3: P_CL - P_SS, scaled by cc
        Dq(1,1) = 1;
        Dq(2,2) = -1;
        if (nest)
            Dq(2,3) = -1;                        % note (2): bare, NOT -cc
            Dq(3,2) =  1;
            Dq(4,1) = -1;  Dq(4,3) = 1;          % note (2): bare, NOT cc
        else
            Dq(2,3) = -cc;
            Dq(3,1) = -1;  Dq(3,2) = 1;
            Dq(4,3) =  cc;
        end
        B = [1;0;-1;0];
        if (has4)                                % resonant vent: CL <-> clvtgt
            vt = 3; if (isfield(pa,'clvtgt')), vt = pa.clvtgt; end
            Df(4,4) =  1;  Df(4,vt) = -1;
            Dq(4,4) =  1;  Dq(vt,4) = -1;
        end
end

C = struct('nch',nch,'ndof',ndof,'Df',Df,'Dq',Dq,'mu',mu,'B',B,'dref',dref);
C.reciprocal = isequal(size(Dq), size(Df')) && norm(Dq - Df', inf) < 1e-12;
end

% NOTE ON FIXING (2). Making nested reciprocal means either applying cc to the
% nested Dq entries, or dropping it from the nested Df pickup. Which is correct
% depends on what clcouple is FOR: it is the m=3b->m=4 reduction knob, and that
% reduction runs non-nested, so nested+clcouple has no established meaning. The
% honest fix is probably to reject clcouple~=1 under nested rather than guess.
% Not done here -- Stage 1 must not change behaviour.
