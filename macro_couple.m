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
if (use3), dref(2) = true; end
% d3 FRAME IS SPLIT BY CHAMBER COUNT, and that is a KNOWN INCONSISTENCY carried
% deliberately (2026-07-29, SN: "split the change: keep the m<4 half, revert
% m=4"). It is exactly the "micromechanics property wearing a macromechanics
% key" this refactor exists to remove, so it should not survive.
%
%   m<4  LAB frame (a3 = s3/m5).  The BM-accelerating frame needs a fictitious
%        force -m5*a1 the force law never supplied. Harmless while d3 was inert;
%        restoring the reaction made it a live unbalanced acceleration and the
%        march diverged unconditionally (k5 x0.001 only delayed it 118 -> 1255;
%        nimp 2->32 changed nothing). Lab frame measured STABLE with d3 live.
%   m>=4 BM frame, unchanged from HEAD, because the same correction makes m=4
%        diverge for a reason not yet identified (see the note in micro26's
%        m>=4 branch). Reverted rather than left broken.
if (numel(dref)>=3), dref(3) = (m>=4); end
% d3 IS IN THE LAB FRAME (2026-07-29, adopted by SN after measurement).
%
% dref(3) used to be true, putting d3 in the BM-ACCELERATING frame
% (a3 = a1 + s3/m5). That frame's equation of motion carries a fictitious force
% -m5*a1, which the force law never supplied, so the frame and the force law
% disagreed. While d3 had no reaction on d1 that cost nothing -- d3 was inert
% and the error had no path to anything. Restoring the reaction (Newton's third
% law) made it a live unbalanced acceleration, and the march diverged
% UNCONDITIONALLY: shrinking k5 1000x only delayed the blowup (sample 118 ->
% 1255), and nimp 2->32 changed nothing (divergence converged to sample 113,
% so the fluid iteration converges -- to an unstable march).
%
% The lab frame is what m5*a3 = -k5*(d3-d1) - r5*(v3-v1) actually says, and it
% measured STABLE with d3 finally live (rel-to-2DOF 9.698 vs 0.000e+00 before).
%
% dref(2) is UNCHANGED and carries the same open question. It is left alone
% because changing it would break the m=1 == m=2 identity, which is a separate
% physics decision rather than a consequence of this one.
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
      if (isfield(pa,'rlsplit') && pa.rlsplit)
        % ---- RL-SPLIT TOPOLOGY (2026-07-29, SN) ---------------------------
        % Chambers are [ST, CL, SV] -- SS IS GONE, replaced by CL. The
        % partition becomes TWO membranes with a fluid space between them:
        %
        %     SV  |  RL (d3)  |  CL  |  BM (d1)  |  ST
        %
        %   d3 (RL) driven by  P_SV - P_CL
        %   d1 (BM) driven by  P_CL - P_ST
        %
        % Sign convention matches the legacy row below, (below - above), so
        % d1's ST term stays +1 exactly as in [1 0 -1].
        %
        % WHY SS IS DROPPED RATHER THAN FROZEN: d2 is SS's ONLY partition
        % coupling (the legacy Df(2,:) = [0 -1 1]). Clamp the TM and SS retains
        % nothing, leaving an UNCOUPLED fluid compartment -- the measured
        % clcouple=0 arm at m=4 diverged at sample 3 exactly that way, the
        % fastest failure of anything tested this session.
        %
        % WHY THIS IS THE PREREQUISITE FOR pa.hbrl: without its own fluid
        % loading the RL is only spring-coupled to the BM through k5/r5, so it
        % FOLLOWS it -- measured max|d3-d1|/max|d1| = 0.1457. A drive
        % proportional to d3 is then nearly IN PHASE with d1, a stiffness change
        % rather than energy injection. micro26's own note says exactly this:
        % "Energy injection is set by the PHASE of the active force against BM
        % velocity". Giving the RL its own pressure difference is what lets it
        % move out of phase.
        %
        % Liu & Neely (2010) is the 1-chamber limit of this: continuity is
        % driven from the RL (eq 9, d_x U = w*xi_r_dot) while pressure acts on
        % the BM (eq 10). That asymmetry is an artifact of collapsing CL; with
        % CL present each membrane carries its own difference and the operator
        % is reciprocal.
        % *** DOES NOT RUN. DEFAULT OFF. *** Diverges at sample 26 for EVERY CL
        % area tested (chsz(2) = 0.05 to 2.00, a 40x span, moved the divergence
        % point by ONE sample). So it is not CL thickness and not a stiffness or
        % CFL limit -- those scale with the parameter. Compare: the reaction
        % instabilities blew up at samples 113-1255 and scaled with coupling
        % strength; an uncoupled chamber (clcouple=0) went at sample 3. A blowup
        % at 26 that is invariant under a large parameter change is the signature
        % of an ILL-POSED system.
        % RULED OUT: CL area; volume conservation (every Dq column sums to zero,
        % as does B); reciprocity (Dq == Df' by construction here).
        % NOT RULED OUT: m=3-specific structure elsewhere that still assumes the
        % old ST/SS/SV roles -- macro26 builds the fluid operator from
        % chsz(1..3) presuming what each chamber adjoins, and there may be
        % boundary or constraint handling keyed to that arrangement beyond
        % Df/Dq. That is where to look next.
        Df(1,:) = [ 1 -1  0];                    % d1 (BM): P_ST - P_CL
        Df(2,:) = [ 0  0  0];                    % d2 (TM): clamped, no coupling
        Df(3,:) = [ 0  1 -1];                    % d3 (RL): P_CL - P_SV
        Dq = Df.';                               % RECIPROCAL BY CONSTRUCTION
        B = [1;0;-1];                            % stapes still drives SV vs ST
      else
        Df(1,:) = [ 1  0 -1];                    % d1: P_ST - P_SV
        Df(2,:) = [ 0 -1  1];                    % d2: P_SV - P_SS
        % d3 (m=3b) is INTERNAL: no fluid compartment, so Df(3,:) stays zero.
        Dq(1,1) =  1;
        Dq(2,2) = -1;
        Dq(3,1) = -1;  Dq(3,2) = 1;
        B = [1;0;-1];
      end
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
