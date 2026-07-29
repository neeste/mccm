function [ss,ii,gam,ohcp,ohcbm] = micro26(pa, cp, st)
%MICRO26  Micromechanics: the partition DOF system.
%
% STAGE 4 of the macro/micro separation (capstone design note: "Micromechanics:
% the partition DOF system, selected by DOF count"). Extracted VERBATIM from
% tdm26/force_cp -- no logic changed, so arch_gate must still read 9/9.
%
% CONTRACT
%   in   pa  model parameters (m, dof, d3int, gam, ohcgain, m3form, ...)
%        cp  place-dependent coefficients from imped26 (k1..r5, m1/m2/m5, hb, gm)
%        st  current state, st.d and st.v over dof*n
%   out  ss  (n x dof) internal partition forces, BEFORE the fluid pressure term
%            (fold_p subtracts Df*p; xpnd_q injects Dq*(mu.*ss))
%        ii  (n x dof) state indices per DOF
%        gam, ohcp, ohcbm  gain and OHC power diagnostics
%
% ss is the interface: everything the macromechanics needs from the partition,
% and the only thing it needs. That is what makes chamber count and DOF count
% independent -- see macro_couple.m for the other side.
%
% KNOWN GAP (design-note decision 1): for m<3 the force law uses an ALGEBRAIC
% shear d3=d1-d2 and defines no k_act/r_act, so it cannot drive a dynamical
% third DOF. Requesting dof>=3 there raises tdm26:noMicro3 rather than failing
% obscurely. That algebraic shear is precisely the variable the design note
% wants promoted to its own mass.

n = pa.n; i1 = 1:n; i2 = (n+1):(2*n); gam = cp.gm; ohcp = 0; ohcbm = 0;
if (pa.m>=4 || (isfield(pa,'d3int') && pa.d3int)), i3=(2*n+1):(3*n); ii=[i1(:) i2(:) i3(:)]; else, ii=[i1(:) i2(:)]; end
if (pa.dof>=4), i4=(3*n+1):(4*n); ii=[ii i4(:)]; end   % vent flow state
if (pa.m<1), ss=0; return; end
d1 = st.d(i1); v1 = st.v(i1); d2 = st.d(i2); v2 = st.v(i2);
if (pa.m>=4 || (isfield(pa,'d3int') && pa.d3int)), dc = st.d(i3); vc = st.v(i3); end   % DOF-3: OC height
if (pa.dof>=4), dq = st.d(i4); vq = st.v(i4); end   % DOF-4: vent flow (integrated)

if (pa.hbnl)
    d3 = cp.hb(:,1) .* d1 + cp.hb(:,2) .* d2;
    if (pa.mmeq == 1), hbt = max(abs(d3) / pa.hbmx,1); gam = cp.gm ./ (1 + pa.hbsc * log(hbt));
    elseif (pa.mmeq == 9), dbt = abs(d3) / pa.hbmx; if (dbt > 1), gam = cp.gm / (1 + pa.hbsc * log(dbt)); end; end
end

if (pa.m<3)
    d3 = d1 - d2; v3 = v1 - v2;
    s1tmp = cp.k1 .* d1 + cp.r1 .* v1; s2tmp = cp.k2 .* d2 + cp.r2 .* v2;
    s3tmp = cp.k3 .* d3 + cp.r3 .* v3; s4tmp = cp.k4 .* d3 + cp.r4 .* v3;
    s1 = -(s1tmp + s3tmp .* cp.gh - s4tmp .* gam); s2 = -(s2tmp - s3tmp);
elseif (pa.m==3)
    % pa.m3form selects the MICROMECHANICS while leaving the 3-chamber
    % HYDRODYNAMICS untouched -- the separation the swap test could not make,
    % because that varied chamber count and parameter set together.
    %   0 (default): FDM-translation form; active force driven by d2 ALONE
    %   1          : m=2-style form; active force driven by the RELATIVE
    %                displacement d3 = d1 - d2
    % Energy injection is set by the PHASE of the active force against BM
    % velocity, so the driving coordinate is precisely what decides whether the
    % force amplifies or dissipates.  Note the "Removes the erroneous k_act*d1
    % stiffening" comment below: this branch was at some point changed FROM a
    % d1-bearing form TO d2-only, which is the very edit under test here.
    % *** pa.m3form IS NOW CONTINUOUS: alpha in [0,1], not a binary switch. ***
    % The two laws differ in TWO places, not one -- the active drive coordinate
    % AND the passive row-2 coupling (-k2/-r2 vs -k3/-r3) -- so a single-
    % coordinate blend cannot reach both endpoints. Interpolating the two zk
    % OPERATORS does, and reduces EXACTLY to each law at alpha = 0 and 1:
    %   s1 active term : k_act*[alpha*d1 + (1-2*alpha)*d2]   (+ same for r_act,v)
    %   s2 row-1 coeff : -[(1-alpha)*k2 + alpha*k3]*d1       (+ same for r,v)
    % alpha=0 -> k_act*d2  and -k2*d1   (FDM d2-only form, the default)
    % alpha=1 -> k_act*(d1-d2) and -k3*d1  (m=2 relative-displacement form)
    % Making the force law CONTINUOUS turns a discrete model-class switch into a
    % fittable parameter -- the reformulation SN's marginal-improvement principle
    % needs, since the chamber count itself cannot be made continuous (chsz->0 is
    % singular in fdm26 and unstable in tdm26).
    % MEASURED TRADE the interpolation is meant to exploit: alpha=0 gives better
    % tuning (maperr 499 vs 1869) while alpha=1 gives a sharper tip (contrast
    % 11.1 vs 9.2) and pushes the apical degeneracy an octave higher.
    al = 0; if (isfield(pa,'m3form')), al = pa.m3form; end
    k_act = cp.gh .* cp.k3 - gam .* cp.k4;
    % ACTIVE RESISTANCE.  fdm26 builds za = z4 = k4/s + r4 (fdm26.m:960) and
    % subtracts zg = gam*za from the row-1 coefficient of V2 (fdm26.m:869,890),
    % so r_act must carry -gam*r4 exactly as k_act carries -gam*k4.  This is a
    % REAL discrepancy vs fdm26 but is currently INERT: pa.r4o = 0 in every
    % parameter set (modpar26c3.m:43, modpar26c4.m:43, par_CEL16.m:36), so
    % cp.r4 == 0 and the term vanishes.  Kept for any parameter set with r4o~=0.
    % (An earlier claim that this omission made the force energy-neutral was
    % WRONG: these are CROSS-coupling terms -- k4*d2 acting on v1 -- which do
    % net work whenever d2 and v1 differ in phase.  m=2 amplifies +39 dB with
    % r4 = 0, which settles it.)
    r4a = 1; if (isfield(pa,'r4act')), r4a = pa.r4act; end
    r_act = cp.gh .* cp.r3 - r4a .* gam .* cp.r4;

    % Interpolated translation of the FDM impedance matrix zk (see alpha note).
    % At alpha=0 this is the exact d2-only form: row 1 carries NO k_act*d1 term
    % (what the older comment called "the erroneous k_act*d1 stiffening"); that
    % term is legitimate only in the alpha>0 relative-displacement law.
    kmix = (1-al) .* cp.k2 + al .* cp.k3;
    rmix = (1-al) .* cp.r2 + al .* cp.r3;

    % Row 1: [z1 + alpha*z_act]*V1 + [(1-2*alpha)*z_act]*V2
    s1 = -(cp.k1 .* d1 + cp.r1 .* v1 ...
           + al .* (k_act .* d1 + r_act .* v1) ...
           + (1-2*al) .* (k_act .* d2 + r_act .* v2));

    % Row 2: -[(1-alpha)*z2 + alpha*zh]*V1 + (z2+zh)*V2
    s2 = -(-kmix .* d1 - rmix .* v1 + (cp.k2 + cp.k3) .* d2 + (cp.r2 + cp.r3) .* v2);
    if (isfield(pa,'d3int') && pa.d3int)
        % m=3b: d3 is an INTERNAL DOF (no fluid compartment). The OHC acts as an
        % internal pair between BM and RL, driven by the shear d2 exactly as in
        % m=3. THE SIGN IS FIXED BY THE REDUCTION GATE: s1 above is m=3's
        % -(k1*d1 + r1*v1 + act), so with d3 frozen the model reproduces m=3
        % exactly; the BM therefore keeps -act and d3 takes the opposite +act.
        % The existing m=4 has this pair INVERTED (+act on BM), which is the
        % leading candidate for its dead amplifier.
        act3 = k_act .* d2 + r_act .* v2;
        f3 = 1; if (isfield(pa,'ohcgain')), f3 = pa.ohcgain; end
        g3 = 1; if (isfield(pa,'ohcsgn')),  g3 = pa.ohcsgn;  end
        act3 = g3 * f3 .* act3;
        s3 = -(cp.k5 .* dc + cp.r5 .* vc - act3);   % RL: +act3 (the reaction)
        bpow = act3 .* v1;  ohcbm = -sum(bpow(isfinite(bpow)));
        wpow = act3 .* (vc - v1); ohcp = -sum(wpow(isfinite(wpow)));
    end
elseif (pa.m>=4)
    % 4-chamber topology, chsz = [ST SS SV CL]. Each DOF is the partition
    % between two chambers, appearing in exactly those two rows with opposite
    % signs (see the qx assignments in cochlea):
    %   d1  basilar membrane   ST <-> SV   also carries the stapes drive
    %   d2  TM-RL shear        SV <-> SS   drives the active force (act)
    %   d3  OC height          SS <-> CL   the OHC acts across this
    % giving a LINEAR CHAIN  ST -d1- SV -d2- SS -d3- CL, with CL a dead-end
    % side compartment (it appears in one row only and sees d3 alone).
    % NOTE d1 was briefly coupled ST<->CL during development; that was reverted
    % but this comment was not updated until 2026-07-25. Verify against the qx
    % assignments, not this header, if they ever disagree again.
    % ANATOMY: the OHC base sits on the Deiters cell (BM side) and its apex is
    % embedded in the RL, so its somatic force acts BETWEEN BM and RL. That is an
    % INTERNAL FORCE PAIR on the OC-height coordinate d3, with the equal-and-
    % opposite reaction on the BM -- NOT a fractional split between two
    % independent forces (one internal force inherently pushes both ends, and the
    % net force on the BM+RL pair is zero, as it must be).
    %   pa.ohcgain scales the pair (1 = full, 0 = passive / no OHC force)
    %   pa.ohcsgn  selects the sense of the pair (+1/-1); which sign AMPLIFIES
    %              rather than damps is settled empirically, not by argument.
    k_act = cp.gh .* cp.k3 - gam .* cp.k4;
    r4a = 1; if (isfield(pa,'r4act')), r4a = pa.r4act; end
    r_act = cp.gh .* cp.r3 - r4a .* gam .* cp.r4;   % ACTIVE RESISTANCE: see the
                             % m==3 branch -- -gam*r4 is what lets the OHC force
                             % do net work (fdm26.m:960 za=k4/s+r4).
    act   = k_act .* d2 + r_act .* v2;      % somatic force from bundle deflection
    fsp = 1;  if (isfield(pa,'ohcgain')), fsp = pa.ohcgain; end
    sgn = 1;  if (isfield(pa,'ohcsgn')),  sgn = pa.ohcsgn;  end
    act = sgn * fsp .* act;
    % ENERGY DIAGNOSTIC. The OHC contributes -act to the height-coordinate force,
    % and (dc-d1) is the relative (BM-RL separation) coordinate, so its power is
    % (-act)*(vc-v1).  ohcp>0 => the force opposes damping and INJECTS energy
    % into the micromechanics; ohcp<0 => it dissipates (wrong sign OR wrong k5
    % phase -- these look identical in WNR magnitude but opposite here).
    wpow = act .* (vc - v1);          % boundary points can give 0*Inf -> NaN
    ohcp = -sum(wpow(isfinite(wpow)));
    % BM WORK. ohcp above is the work on the RELATIVE (BM-to-RL) coordinate,
    % which is the right measure for the internal force pair. It does NOT say
    % whether the BM is being amplified. The BM receives +act (the reaction in
    % s1 below), so the power delivered to BM motion is act.*v1. ohcbm>0 means
    % the OHC force opposes BM damping and drives the travelling wave; ohcbm~0
    % with ohcp>0 means the force is injecting energy into the OC-height
    % coordinate only and never reaching the BM, which would explain an
    % amplifier that stays near 2.5 dB regardless of gain or damping.
    bpow  = act .* v1;
    ohcbm = sum(bpow(isfinite(bpow)));
    % SIGN FIXED BY THE m=3b REDUCTION GATE. m=3b has s1 = -(k1*d1 + r1*v1 +
    % act) so the BM receives -act, and s3 = -(k5*dc + r5*vc - act) so d3
    % receives +act. m=4 must reduce to that at clcouple=0, so it carries the
    % SAME pair. The legacy m=4 had it INVERTED (+act on the BM); the bisection
    % confirmed that was the cause of its dead amplifier (correcting the sign
    % moved gain from +3.43 to ~+78 dB). pa.m4legacy=1 restores the old form.
    lg = 0; if (isfield(pa,'m4legacy')), lg = pa.m4legacy; end
    if (lg)
        s1 = -(cp.k1 .* d1 + cp.r1 .* v1 - act);      % BM   : legacy (inverted)
        s3 = -(cp.k5 .* dc + cp.r5 .* vc + act);      % OC height : legacy
    else
        s1 = -(cp.k1 .* d1 + cp.r1 .* v1 + act);      % BM   : -act (matches m=3b)
        s3 = -(cp.k5 .* dc + cp.r5 .* vc - act);      % OC height : +act reaction
    end
    s2 = -(-cp.k2 .* d1 - cp.r2 .* v1 + (cp.k2 + cp.k3) .* d2 + (cp.r2 + cp.r3) .* v2);
end
s1(n) = 0; s2(n) = 0;
% THIRD DOF. The old one-liner did `s3(n)=0` unconditionally, which AUTO-CREATES
% s3 as a 1xn ROW when no force law has defined it -- then [s1 s2 s3] fails with
% "Dimensions of arrays being concatenated are not consistent", several frames
% from the actual cause. That is what nch=1,dof=3 hit. Diagnose it properly.
need3 = (pa.m>=4) || (isfield(pa,'d3int') && pa.d3int);
if (need3)
    if (~exist('s3','var') || numel(s3) ~= n)
        error('tdm26:noMicro3', ...
            ['A third DOF was requested (m=%d, d3int=%d) but no third-DOF force ' ...
             'law ran.\nThe d3int block lives inside the m==3 branch. For m<3 the ' ...
             'force law uses an ALGEBRAIC shear d3=d1-d2 (line 797) and defines ' ...
             'no k_act/r_act, so there is nothing to drive a dynamical third DOF.' ...
             '\nThis is design-note DECISION 1, not a plumbing gap: is DOF3 the ' ...
             'OHC/cilia site promoted from that algebraic shear to its own mass, ' ...
             'or a different partition of the three? That choice defines the m<3 ' ...
             'third-DOF force law and cannot be guessed by a refactor.'], ...
            pa.m, double(isfield(pa,'d3int') && pa.d3int));
    end
    s3 = s3(:); s3(n) = 0; ss = [s1 s2 s3];
else
    ss = [s1 s2];
end
if (pa.dof>=4)
    % VENT restoring force. The pressure drive (P_CL - P_target) is added in
    % fold_p, exactly as d1/d2/d3 receive theirs, so s4 here carries ONLY the
    % channel's own stiffness and resistance. With clvk = clvr = 0 this is
    % identically zero and the vent reverts to the pure inertance that the a2
    % stamp already represents -- that is the reduction gate.
    s4 = -(cp.clvk .* dq + cp.clvr .* vq);
    s4(n) = 0; ss = [ss s4];
end
end

%==========================================================
% VISUALIZATION & PROTOCOLS
%==========================================================

