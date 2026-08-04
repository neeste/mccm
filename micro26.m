function [ss,ii,gam,ohcp,ohcbm,adrv] = micro26(pa, cp, st)
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
% shear dh=d1-d2 and defines no k_act/r_act, so it cannot drive a dynamical
% third DOF. Requesting dof>=3 there raises tdm26:noMicro3 rather than failing
% obscurely. That algebraic shear is precisely the variable the design note
% wants promoted to its own mass.

n = pa.n; i1 = 1:n; i2 = (n+1):(2*n); gam = cp.gm; ohcp = 0; ohcbm = 0;
adrv = zeros(n,1);   % instantaneous OHC active drive, for the RC pole (pa.ohctau)
% THE RC POLE IS IMPLEMENTED FOR THE m>=3 FORCE LAW ONLY. The m<3 branch carries
% a different active term -- cp.gh.*gam.*s4tmp, with gh multiplying the whole
% active force rather than sitting inside k_act -- so filtering it correctly
% needs its own drive definition. Refusing is better than silently applying the
% pole to one branch and not the other, which would make m=1 and m=3 describe
% different physics under the same flag.
if (isfield(pa,'ohctau') && ~isempty(pa.ohctau) && pa.ohctau > 0 && pa.m < 3)
    error('micro26:ohctau_m', ...
          'pa.ohctau (OHC RC pole) is implemented for m>=3 only; got m=%d.', pa.m);
end
% THIRD-DOF PREDICATE, defined ONCE. This condition previously appeared at three
% separate sites (the i3 index, the dc/vc extraction, and the ss assembly) and
% adding the m<3 case to only one of them produced "Unrecognized function or
% variable 'dc'" -- the index and the force law disagreed about whether a third
% DOF existed. Duplicated predicates are how that happens; keep this single.
has3 = (pa.m>=4) || (isfield(pa,'d3int') && pa.d3int) || ...
       (pa.m<3 && isfield(pa,'dof') && pa.dof>=3);
if (has3), i3=(2*n+1):(3*n); ii=[i1(:) i2(:) i3(:)]; else, ii=[i1(:) i2(:)]; end
if (pa.dof>=4), i4=(3*n+1):(4*n); ii=[ii i4(:)]; end   % vent flow state
if (pa.m<1), ss=0; return; end
d1 = st.d(i1); v1 = st.v(i1); d2 = st.d(i2); v2 = st.v(i2);
if (has3), dc = st.d(i3); vc = st.v(i3); end   % DOF-3: OC height / cortilymph
if (pa.dof>=4), dq = st.d(i4); vq = st.v(i4); end   % DOF-4: vent flow (integrated)

% BUNDLE COORDINATE, DEFINED ONCE for every branch below (2026-07-29).
% dhb/vhb are the hair bundle xi_h and its velocity, via cp.hb -- see hbmix.m.
% Defined here rather than inside each branch because the has3 predicate taught
% what duplicated definitions cost: it lived at four sites, one was updated, and
% the result was "Unrecognized function or variable 'dc'".
%
% THE ACTIVE DRIVE IS NOW THE BUNDLE, not d2. At m>=3 the default hb = [0 1 0]
% makes dhb == d2 exactly, so this is BIT-IDENTICAL by default at every chamber
% count. It differs ONLY under pa.hbrl, where hb = [0 -1 gh] gives
% dhb = gh*d3 - d2; with the TM clamped (pa.tmrigid) that becomes gh*d3, the
% RL-driven MET of Liu & Neely (2010), i_r = I(alpha_v*xi_r_dot + alpha_d*xi_r).
%
% WHY THIS WAS BLOCKING: with the drive on d2, clamping the TM sets act == 0 and
% the amplifier vanishes -- measured, max|d1| fell from 8.03e-03 to 1.13e-08.
% The clamp is only usable once the drive no longer depends on the clamped DOF.
if (has3), dhb = hbmix(cp.hb, d1, d2, dc); vhb = hbmix(cp.hb, v1, v2, vc);
else,      dhb = hbmix(cp.hb, d1, d2);     vhb = hbmix(cp.hb, v1, v2);
end

if (pa.hbnl)
    % dh IS THE HAIR BUNDLE (xi_h), an ALGEBRAIC shear, renamed from d3 on
    % 2026-07-29. It was called d3 while the DYNAMICAL third DOF was also called
    % d3 in comments (its state is dc/vc), and that collision cost real
    % confusion. dh matches cochlea_proc.docx's xi_h. See macro26.m:31 for what
    % the shear IS: gh*d1 - d2 at m<3 (d2 = TM), and just d2 at m>=3 (d2 = the
    % bundle itself, TM eliminated).
    dh = dhb;   % same bundle coordinate; computed once above
    if (pa.mmeq == 1), hbt = max(abs(dh) / pa.hbmx,1); gam = cp.gm ./ (1 + pa.hbsc * log(hbt));
    elseif (pa.mmeq == 9), dbt = abs(dh) / pa.hbmx; if (dbt > 1), gam = cp.gm / (1 + pa.hbsc * log(dbt)); end; end
end

% MICROMECHANICS VARIANT (2026-07-28, SN): "when m<3 the micro-mechanics should
% mimic m=4 via interaction between SS and CL without longitudinal coupling."
% At m3form=0 the m=3 force law is ALREADY identical to m=4's -- both give
% s1 = -(k1*d1 + r1*v1 + k_act*d2 + r_act*v2) and the same shear row -- and
% d3int is m=4's d3 mechanism with NO fluid compartment, i.e. no longitudinal
% coupling. So "mimic m=4" reduces to: use the m=3 branch when a third DOF is
% requested, whatever the chamber count. This is the FIRST behaviour change of
% the refactor; m=1/m=2 with dof=2 are untouched and still take the legacy law.
use3 = (pa.m==3) || (pa.m<3 && isfield(pa,'dof') && pa.dof>=3);
if (pa.m<3 && ~use3)
    % pa.ghlever: gh IS A LEVER (SN, 2026-07-29), so it belongs on BOTH the
    % displacement and the force. The legacy law applies it inconsistently --
    % THREE different treatments of gh in this one branch:
    %   cp.hb      gh*d1 - d2      lever on DISPLACEMENT   (used by the MET gain
    %                              and, since today, the active drive)
    %   passive    d1 - d2, with gh on the FORCE (s3tmp.*cp.gh)
    %   active     d1 - d2, with NO gh at all (-s4tmp.*gam)
    % A lever ratio means displacement gh*d1 at the bundle and force gh*F back at
    % the BM, so the correct law carries gh in both places for BOTH terms. That
    % makes the force law agree with cp.hb, which has had it right all along.
    %
    % ADOPTED UNCONDITIONALLY, and it is FREE TODAY: cp.gh is EXACTLY 1 at m<3
    % (measured, max|gh-1| = 0), so m=1 and m=2 are bit-identical -- max|d1|
    % 1.9757e-07, maperr 104.629, amp +39.11, unchanged.
    %
    % But gh is NOT 1 at m>=3 (max|gh-1| = 4.4e-03, from gpo*exp(gpe*x+gpq*q)).
    % So the inconsistency was invisible ONLY because this branch never ran where
    % gh departs from unity. Substituting this law into m=3 -- the stated next
    % step -- makes gh vary, at which point the three treatments give three
    % different force laws. Fixing it now costs nothing and removes that trap.
    dh = cp.gh .* d1 - d2;  vh = cp.gh .* v1 - v2;     % lever on DISPLACEMENT
    s1tmp = cp.k1 .* d1 + cp.r1 .* v1; s2tmp = cp.k2 .* d2 + cp.r2 .* v2;
    s3tmp = cp.k3 .* dh + cp.r3 .* vh; s4tmp = cp.k4 .* dh + cp.r4 .* vh;
    s1 = -(s1tmp + cp.gh .* (s3tmp - gam .* s4tmp));   % lever on FORCE, both terms
    s2 = -(s2tmp - s3tmp);   % TM sits at the bundle end: no lever here
elseif (use3)
    % pa.m3form selects the MICROMECHANICS while leaving the 3-chamber
    % HYDRODYNAMICS untouched -- the separation the swap test could not make,
    % because that varied chamber count and parameter set together.
    %   0 (default): FDM-translation form; active force driven by d2 ALONE
    %   1          : m=2-style form; active force driven by the RELATIVE
    %                displacement dh = d1 - d2
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
    % m3form SPLIT INTO TWO INDEPENDENT CHOICES (2026-07-29, SN). One alpha used
    % to control two unrelated things, which is why "substitute the TM law" and
    % "drive from the bundle" pulled opposite directions on one knob:
    %
    %   m3row  WHICH COORDINATE IS SLAVED -- the s2 row coupling.
    %          0 -> kmix = k2   d2 behaves as the BUNDLE   (current m>=3)
    %          1 -> kmix = k3   d2 behaves as the TM       (the m<3 law)
    %
    %   m3drv  WHICH COORDINATE DRIVES the active force.
    %          0 -> k_act*dhb          the bundle ITSELF   (what pa.hbrl needs)
    %          1 -> k_act*(d1 - dhb)   a difference from d1 (the m<3 form)
    %
    % SN's framing makes these orthogonal: the two 2-DOF models are the two ways
    % of slaving one coordinate of a 3-DOF system (d3 = gh*d1, or d2 = eps*d3),
    % while the drive coordinate is a separate question about where MET reads.
    % The 2010 target is m3row=1 (d2 as TM) WITH m3drv=0 (bundle drive) -- a
    % combination the single knob could not express.
    %
    % BACKWARD COMPATIBLE: both default to m3form, so setting m3form alone
    % reproduces the old coupled behaviour exactly, including m3form=1 being the
    % m<3 law to roundoff.
    al = 0; if (isfield(pa,'m3form')), al = pa.m3form; end
    ar = al; if (isfield(pa,'m3row')), ar = pa.m3row; end   % slaving / s2 row
    ad = al; if (isfield(pa,'m3drv')), ad = pa.m3drv; end   % drive coordinate
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
    kmix = (1-ar) .* cp.k2 + ar .* cp.k3;
    rmix = (1-ar) .* cp.r2 + ar .* cp.r3;

    % ---- OHC RC POLE (pa.ohctau > 0; 0 = off = the legacy instantaneous law) ----
    % Za = Ka/(i*omega) is PURE STIFFNESS (which is why r4o = 0 in every
    % parameter set -- the formulation, not an omission), so the active element
    % is negative STIFFNESS: 90 deg. One RC pole on the OHC lateral-wall voltage
    % adds a second 90 deg above its corner, and 180 deg turns the element into
    % frequency-dependent negative DAMPING -- the cochlear-amplifier mechanism
    % (SN's tutorial, cochamp.htm third scene: below ~1 kHz the voltage follows
    % the current and the contraction merely absorbs BM motion; above it the
    % contraction lags 90 deg and pumps energy in, "rocking a boat").
    %
    % cur.vohc / cur.dvohc have been allocated since tdm24 and never read or
    % written; this is what they are for. vohc holds the FILTERED active drive,
    % so it carries the same units as adrv below and the gain gam still
    % multiplies it. tdm26 advances it ONCE per step -- see ohc_rc_step there,
    % and the note explaining why it cannot be advanced inside accel.
    %
    % fdm26 gets the SAME pole as z4 -> z4/(1+s*tau) inside its own imped, which
    % is the single place both the 1-/2-chamber Yb and the multi-chamber zg read
    % za from. Both solvers must carry it or maperr silently compares two
    % different models.
    uact = ad .* d1 + (1-2*ad) .* dhb;      % active drive coordinate
    wact = ad .* v1 + (1-2*ad) .* vhb;
    adrv = cp.k4 .* uact + r4a .* cp.r4 .* wact;   % instantaneous drive, pre-gain
    rcon = isfield(pa,'ohctau') && ~isempty(pa.ohctau) && pa.ohctau > 0;

    % Row 1: [z1 + alpha*z_act]*V1 + [(1-2*alpha)*z_act]*V2
    if (rcon && isfield(st,'vohc') && ~isempty(st.vohc))
        % Passive part instantaneous; ACTIVE part is the filtered OHC voltage.
        % Algebraically identical to the legacy line when vohc == adrv:
        %   k_act*uact + r_act*wact == gh*k3*uact + gh*r3*wact - gam*adrv
        s1 = -(cp.k1 .* d1 + cp.r1 .* v1 ...
               + cp.gh .* cp.k3 .* uact + cp.gh .* cp.r3 .* wact ...
               - gam .* st.vohc);
    else
        % LEGACY EXPRESSION KEPT VERBATIM. The identity above holds in exact
        % arithmetic but NOT bit-for-bit: k_act = gh*k3 - gam*k4 rounds once
        % before multiplying by uact, while the split form rounds twice. Reusing
        % the new expression here would shift the last bits of every result the
        % project has banked, so the off path stays byte-identical by keeping
        % the original line rather than by trusting the algebra.
        s1 = -(cp.k1 .* d1 + cp.r1 .* v1 ...
               + ad .* (k_act .* d1 + r_act .* v1) ...
               + (1-2*ad) .* (k_act .* dhb + r_act .* vhb));
    end

    % Row 2: -[(1-alpha)*z2 + alpha*zh]*V1 + (z2+zh)*V2
    s2 = -(-kmix .* d1 - rmix .* v1 + (cp.k2 + cp.k3) .* d2 + (cp.r2 + cp.r3) .* v2);
    % d3int, or an explicit third DOF at m<3 (which implies the internal form:
    % there is no CL chamber below m=4 for it to couple to).
    if (has3)
        % m=3b: d3 is an INTERNAL DOF (no fluid compartment). The OHC acts as an
        % internal pair between BM and RL, driven by the shear d2 exactly as in
        % m=3. THE SIGN IS FIXED BY THE REDUCTION GATE: s1 above is m=3's
        % -(k1*d1 + r1*v1 + act), so with d3 frozen the model reproduces m=3
        % exactly; the BM therefore keeps -act and d3 takes the opposite +act.
        % The existing m=4 has this pair INVERTED (+act on BM), which is the
        % leading candidate for its dead amplifier.
        act3 = k_act .* dhb + r_act .* vhb;
        f3 = 1; if (isfield(pa,'ohcgain')), f3 = pa.ohcgain; end
        g3 = 1; if (isfield(pa,'ohcsgn')),  g3 = pa.ohcsgn;  end
        act3 = g3 * f3 .* act3;
        % PASSIVE ATTACHMENT IS RELATIVE (2026-07-28, SN). k5/r5 connect the
        % OHC/OC-height DOF to the BM, so they act on (d3 - d1), not on d3
        % alone. The old form used the ABSOLUTE dc, which is a spring to
        % GROUND, while the acceleration was already referenced to the BM
        % (a(i3) = a(i1) + s3/m5, i.e. m5 acting on the RELATIVE
        % acceleration). Those two were inconsistent: a grounded stiffness
        % with relative inertia.
        % An attachment is a FORCE PAIR, so d1 receives the reaction. That is
        % what gives d3 a path back to the BM -- without it d3 was a driven
        % observer with no back-action, which is why k5/r5/m5 measured as
        % EXACTLY inert in the distillation objective.
        % EXTENSION IS (d3 - d1), and d3 is in the LAB frame (macro_couple sets
        % dref(3)=false for m<4). The two go together and were MEASURED:
        %   lab frame + k5*dc        -> d3 has no drive at all, stays at zero
        %   lab frame + k5*(dc-d1)   -> STABLE, and d3 finally live (rel 9.698)
        % The "d3 is the OC height, so the spring acts on dc itself" reading needs
        % the BM-frame acceleration PLUS a fictitious force -m5*a1 that the force
        % law does not supply, so it is not self-consistent here.
        dcr = dc - d1; vcr = vc - v1;         % BM-to-OHC extension
        spas = cp.k5 .* dcr + cp.r5 .* vcr;   % passive attachment force
        % ACTIVE FORCE REACHES THE BM ONLY WHEN m<4 (SN, 2026-07-28, and
        % cochlea_proc.docx: "an active force gamma*Za*xi_h ... that acts only on
        % BM and is proportional to HB displacement"). It is already in s1 above,
        % via the m=3 row-1 term (1-2al)*(k_act*d2 + r_act*v2) -- equation 24's
        % gh*Zh - gamma*Za. So d3 gets NO active term here; it is driven purely by
        % its passive attachment to the BM, which is what "d1 drives both d2 and
        % d3" means at m<4.
        % At m=4 the force instead acts BETWEEN BM and RL, and the m>=4 branch
        % below keeps that pair. That is the one place the two differ.
        s3 = -spas;                        % OHC DOF: passive attachment only
        s1 = s1 + spas;                    % BM: reaction of that attachment
        % NOTE ohcgain/ohcsgn now have NO effect at m<4: they scaled act3, which
        % no longer reaches d3. They remain live at m=4, where the pair exists.
        % act3 is NOT zeroed -- the diagnostics below still want it. ohcbm
        % (act3.*v1) remains meaningful: the active force does act on the BM, via
        % s1. ohcp (act3.*(vc-v1)) is now MOOT at m<4, since the active force no
        % longer acts across the BM-to-OHC coordinate; read it only at m=4.
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
    act   = k_act .* dhb + r_act .* vhb;    % somatic force from bundle deflection
                                            % (dhb == d2 at default hb, so m=4 is
                                            % unchanged; the comment finally matches)
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
        % pa.m4legacy deliberately keeps the OLD grounded passive element AND the
        % inverted sign pair, so it remains a faithful pre-correction reference.
    else
        % m=4 IS UNCHANGED FROM HEAD, DELIBERATELY (2026-07-29, SN: "split the
        % change: keep the m<4 half, revert m=4").
        %
        % The m<4 correction below/above -- lab frame + reaction -- was applied
        % here too and made m=4 DIVERGE. That was bisected, not guessed: m=4 runs
        % clean at HEAD (max|d1| 1.04e-05), and the reaction alone is responsible.
        % Eliminated by measurement, each in both directions:
        %   frame        diverges in BM frame AND lab frame
        %   extension    diverges with k5*dc AND k5*(dc-d1)
        %   topology     diverges with nested=0, clcouple=0, vent removed
        %   reciprocity  Dq == Df' holds in every configuration (not the cause)
        %   somatic loop removing act from d3 still diverges (1.9e+07)
        %
        % With act off d3, m=4's third-DOF law is IDENTICAL to m=3b's, and m=3b is
        % stable. The one difference left is that m=4's d3 is FLUID-COUPLED. So the
        % reaction is incompatible with a d3 that also exchanges with a chamber,
        % as currently formulated.
        %
        % HYPOTHESIS for whoever picks this up, NOT a finding: xpnd_q injects each
        % DOF's force with mu = [1, m1/m2, m1/m5], weights derived when d3 reached
        % the fluid only through its own chamber. The reaction sends part of that
        % same force into d1, which injects at weight 1, so one force can enter the
        % fluid twice under two weights. m=3b never triggers it (no chamber). That
        % needs mu re-derived for a d3 that is both fluid-coupled and attached,
        % which is a derivation, not a debugging step.
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
need3 = has3;
if (need3)
    if (~exist('s3','var') || numel(s3) ~= n)
        error('tdm26:noMicro3', ...
            ['A third DOF was requested (m=%d, d3int=%d) but no third-DOF force ' ...
             'law ran.\nThe d3int block lives inside the m==3 branch. For m<3 the ' ...
             'force law uses an ALGEBRAIC shear dh=d1-d2 and defines ' ...
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

