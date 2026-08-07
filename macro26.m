function [cp, me] = macro26(pa, me)
%MACRO26  Macromechanics: fluid chambers and their coupling to the partition.
%
% STAGE 4 of the macro/micro separation (capstone design note: "Macromechanics:
% fluid and chambers, selected by chamber count nch"). Extracted VERBATIM from
% tdm26/cochlea -- no logic changed, so arch_gate must still read 9/9.
%
% CONTRACT
%   in   pa  model parameters (m/nch, chsz, nested, clvent, clcouple, ...)
%        me  middle-ear state
%   out  cp  place-dependent coefficients PLUS
%            cp.aa  the block-tridiagonal fluid operator
%            cp.mc  the macro/micro interface from macro_couple (Df, Dq, mu, B, dref)
%
% cp.mc is the whole interface to the micromechanics. Everything else here is
% chambers: areas, the fluid stencil, boundary conditions, chamber-to-chamber
% couplings such as the CL vent.
%
% Calls imped26 for the place-dependent impedance parameters, which belong to
% the MICROmechanics but are built once alongside the chamber geometry.

m = pa.m; n = pa.n;
if (~isfield(pa,'aflom_fac')), pa.aflom_fac = 1; end % default: current tdm26 coupling
cp = imped26(pa);
% Carried on cp because fold_p (which applies it) receives cp, not pa.
% Set here, OUTSIDE the chamber-count branches, so a rigid-BM run is
% honoured for EVERY m. Setting it inside the m>=4 branch would leave
% m=1/2/3 silently running FREE and make the comparison meaningless.
cp.bmrigid = 0; if (isfield(pa,'bmrigid')), cp.bmrigid = pa.bmrigid; end
% pa.tmrigid: KINEMATIC CLAMP of d2 to the LAB frame. Toward Liu & Neely (2010),
% which does not represent the TM at all ("tectorial-membrane motion is not
% considered explicitly"). Clamping is the correct way to express d2 -> 0:
% stiffening k2 instead would create an unresolvable mode, exactly as the
% bmrigid note above records for k1o. Set here, outside the chamber branches, so
% it is honoured for every m.
cp.tmrigid = 0; if (isfield(pa,'tmrigid')), cp.tmrigid = pa.tmrigid; end
cp.clcouple = 1; if (isfield(pa,'clcouple')), cp.clcouple = pa.clcouple; end
% cp.hb IS THREE COLUMNS as of 2026-07-29, over [d1 d2 d3]. The third column is
% ZERO in both legacy forms, so this is behaviour-preserving. Always read it via
% hbmix() -- there were five consumers and a 2-column read is silently wrong
% whenever column 3 is nonzero.
if (m<3)
    cp.hb=[cp.gh -ones(n,1) zeros(n,1)];   % bundle = gh*d1 - d2   (d2 = TM)
else
    cp.hb=[zeros(n,1) ones(n,1) zeros(n,1)]; % bundle = d2         (d2 IS bundle)
end
% pa.hbrl: BUNDLE REFERENCED TO THE RL (d3) instead of the BM (d1). This is Liu &
% Neely (2010)'s micromechanics, where MET current is driven by RL motion
% (i_r = I(alpha_v*xi_r_dot + alpha_d*xi_r)) and the TM is not represented at
% all. With the TM clamped, bundle = gh*d3, driven purely by the RL.
% NOTE this changes only WHERE THE BUNDLE IS MEASURED. The active force's drive
% coordinate in the m>=3 force law is still d2 directly (micro26's
% k_act*d2 + r_act*v2), NOT routed through hb, so referencing the bundle to d3
% does NOT by itself move the active drive. That is a separate change.
if (isfield(pa,'hbrl') && pa.hbrl)
    % COEFFICIENT ON d3 IS 1, NOT gh (corrected 2026-07-29, SN's framing).
    % The bundle is d3 - d2. gh belongs to the SLAVING RELATION d3 = gh*d1, not
    % to the bundle definition: it appears in the m<3 row [gh -1 0] only because
    % substituting gh*d1 for d3 carries it along. Once d3 is an explicit
    % coordinate the lever is no longer implicit, so keeping gh here would
    % DOUBLE-COUNT it -- easy to miss, since gh also multiplies k3 inside k_act.
    %
    % The two 2-DOF models are the two ways to slave one coordinate:
    %   d3 = gh*d1   RL slaved to BM  -> bundle = gh*d1 - d2   (the m<3 row)
    %   d2 = eps*d3  TM slaved to RL  -> bundle = (1-eps)*d3   (2010, eps=0)
    % Neither is the current m>=3 row [0 1 0], which instead REDEFINES d2 to be
    % the bundle. That change of variable, not a slaving, is why d2 means
    % different things above and below m=3.
    % SIGN IS SELECTABLE (2026-08-06). pa.hbrl = +1 keeps bundle = d3 - d2 exactly
    % as before (bit-identical); pa.hbrl = -1 gives d2 - d3. Added because under
    % pa.rlsplit the RL becomes the FLUID-COUPLED membrane and the TM the internal
    % one, which is the reverse of the arrangement this sign was chosen for, and
    % the net OHC work on the BM comes out NEGATIVE there (-7.5e-05 against
    % +2.3e-04 with hbrl off). A mechanism that absorbs where it should inject is
    % the RC-pole failure repeating, and the sign is the first thing to rule out.
    sg = sign(pa.hbrl);
    cp.hb=[zeros(n,1) -sg*ones(n,1) sg*ones(n,1)];   % bundle = sg*(d3 - d2)
end
dx = pa.xl / (n - 1);
if (m < 1), cp.abmom=0; cp.alfx = me.mst / (2 * pa.rho * dx); return; end
cp.abmom = (cp.bw * dx) ./ cp.m1;
cp.alfx = (pa.ast / pa.mst) / cp.abmom(1);
% macro_couple is built HERE, before the fluid operator, because the operator's
% chamber-coupling block is now DERIVED from it rather than written out. All it
% reads (m1/m2/m5/clvm/clvk) comes from imped26 above, so this is safe; the m<1
% early return above still exits without cp.mc, exactly as before.
cp.mc = macro_couple(pa, cp);
aflom = cp.ac / (pa.aflom_fac * pa.rho * dx);
% CHAMBER-SIZE NORMALIZATION. The legacy form rescaled chsz to a FIXED SUM of
% 2, which meant ADDING A CHAMBER RESIZED THE EXISTING ONES: m=4's raw
% [0.95 0.05 1.0 0.05] sums to 2.05, so ST and SV were shrunk 2.44%% before
% anything else happened. That made an EXACT reduction across chamber counts
% IMPOSSIBLE (for m=4's first three to match m=3's, the fourth must be zero,
% which leaves the CL row singular), and near-critical hypersensitivity turned
% that 2.44%% into a 54%% displacement difference in the m=4 reduction gate.
% Now chsz is used AS GIVEN. m=1, m=2 and m=3 are UNAFFECTED because their
% chsz already sums to exactly 2; only m=4 changes. pa.chsznorm=1 restores the
% legacy behaviour. NOTE this makes the OVERALL chsz SCALE a real parameter
% again (it scales the fluid inertance against the partition masses), which
% the fixed-sum constraint had been suppressing.
if (isfield(pa,'chsznorm') && pa.chsznorm)
    pa.chsz = pa.chsz * (2 / sum(pa.chsz));   % legacy fixed-sum form
end
kk = 2:(n-1); mm = m*m;
a1=zeros(n,mm);a2=zeros(n,mm);a3=zeros(n,mm);

if (m<3)
    a3(1,1) = -aflom(1) ./ cp.abmom(1); a2(1,1) =  1 + cp.alfx - a3(1);
    a1(kk,1) = -aflom(kk-1) ./ cp.abmom(kk); a3(kk,1) = -aflom(kk) ./ cp.abmom(kk);
    a2(kk,1) = 1 - a1(kk) - a3(kk); a1(n,1) = -1; a2(n,1) = 1;
elseif m==3
    % CHAMBER COUPLING IS NOW DERIVED, not written out (2026-07-29):
    %     a2(k) = diag(L_p + L_c) + Dq*diag(mu(k,:))*Df
    % macro_shadow verifies the identity (2.842e-14 at m=3 AND m=4). The old
    % hand-written form -- off-diagonals -1 for ST<->SV via d1, -mu for SS<->SV
    % via d2 -- was a SECOND copy of the topology macro_couple already held, and
    % the two could disagree. pa.rlsplit changed Df/Dq and this block did not
    % follow, over-determining the system: divergence at sample 26, INVARIANT
    % under a 40x change in CL area. Deriving it means any topology expressed in
    % Df/Dq is automatically consistent here.
    % mu now enters by DOF INDEX, so each chamber pair carries whichever DOF
    % actually mediates it. The old form froze m1/m2 into chamber 2 because SS
    % happened to be mediated by d2.
    C3 = cp.mc; nc3 = C3.nch;
    dgi = ((1:nc3)-1)*nc3 + (1:nc3);          % column-major diagonal indices
    for k=2:(n-1)
        Lp = aflom(k-1)*pa.chsz(1:nc3)./cp.abmom(k);
        Lc = aflom(k)  *pa.chsz(1:nc3)./cp.abmom(k);
        A2 = diag(Lp(:).' + Lc(:).') + C3.Dq * diag(C3.mu(k,:)) * C3.Df;
        a2(k,:)   = A2(:).';
        a1(k,dgi) = -Lp(:).';
        a3(k,dgi) = -Lc(:).';
    end
    % Basal Boundary Condition (Restores absolute ground)
    L1_c = aflom(1)*pa.chsz(1)./cp.abmom(1);
    L2_c = aflom(1)*pa.chsz(2)./cp.abmom(1);
    L3_c = aflom(1)*pa.chsz(3)./cp.abmom(1);
    a3(1,1) = -L1_c; a3(1,5) = -L2_c; a3(1,9) = -L3_c;        
    % Apply 2*cp.alfx to compensate for ST and SV being in series.
    % The off-diagonals remain -1 to preserve the ambient pressure ground.
    bcx=1; if (isfield(pa,'bcx')), bcx=pa.bcx; end  % DIAGNOSTIC: scale stapes boundary admittance only (drive @L174 untouched); matched ~ 2*alfx=L1_c
    % Same derived block, plus the stapes admittance on the chambers facing the
    % stapes and round window (first and last). That placement is unchanged and
    % remains correct under rlsplit, where the chambers are [ST, CL, SV].
    Lc1 = aflom(1)*pa.chsz(1:nc3)./cp.abmom(1);
    A2b = diag(Lc1(:).') + C3.Dq * diag(C3.mu(1,:)) * C3.Df;
    sa = zeros(1,nc3); sa(1) = 2*cp.alfx*bcx; sa(nc3) = sa(nc3) + 2*cp.alfx*bcx;
    a2(1,:) = reshape(A2b + diag(sa), 1, []);
    % Apical Boundary Condition
    a2(n,:) = 0;
    a1(n,1) = -1; a2(n,1) = 1;
    a1(n,5) = -1; a2(n,5) = 1;
    a1(n,9) = -1; a2(n,9) = 1;
elseif (m>=4)
    % 4 chambers: 1=ST 2=SS 3=SV 4=CL.  Three partition DOFs:
    %   d1 (BM)          couples ST<->SV   coefficient 1     (main drive path)
    %   d2 (TM-RL shear) couples SV<->SS   coefficient mu2 = m1/m2
    %   d3 (OC height)   couples CL<->SS   coefficient mu3 = m1/m5
    % CL is a SIDE compartment pumped by OHC somatic motility via d3; it is NOT
    % in the stapes drive path. Routing d1 through CL (an earlier attempt) made
    % CL a dead end whose only outlet was the series chain d3->SS->d2->SV, which
    % loaded the BM and suppressed d1 ~8x at every CL size, so the compression
    % never engaged and the model ran exactly linear.
    % Each DOF adds +c to both self-terms and -c to the two cross-terms (the same
    % symmetric pattern as the 3-chamber), keeping the operator symmetric and
    % reciprocal with the xpnd_q / fold_p D/D^T pair.
    mu2 = cp.m1 ./ max(cp.m2, 1e-12);
    mu3 = cp.m1 ./ max(cp.m5, 1e-12);
    % pa.clcouple in [0,1] scales d3's FLUID coupling, giving the m=3b -> m=4
    % reduction a CONTINUOUS gate. At 0, d3 has no fluid coupling and CL is a
    % decoupled chamber, so the model must reproduce m=3b EXACTLY; at 1 it is
    % the full 4-chamber. chsz(4)->0 does NOT serve as this limit: there mu3
    % cancels algebraically and d3 still carries a pressure constraint, whereas
    % in m=3b it has none.
    cc = 1; if (isfield(pa,'clcouple')), cc = pa.clcouple; end
    nest = 0; if (isfield(pa,'nested')),  nest = pa.nested; end
    vent = 0; if (isfield(pa,'clvent')),  vent = pa.clvent; end
    % VENT TARGET, as a chamber index: 3 = SV (default), 1 = ST (legacy).
    % SN's decision that CL is CARVED FROM SV makes SV the parent, and a space
    % carved from a pool communicates with the pool it came from. The original
    % ST vent was worse than inert: it stamps G onto the SAME (ST,CL) entries
    % that already carry d1, so it stiffened the BM's own pressure difference
    % rather than opening a second path, which is why raising it drove d1 down.
    % Venting to SV instead restores an ST<->SV exchange route THROUGH CL while
    % leaving p_ST - p_CL intact as the BM drive. Large clvent should then
    % recover the appended amplifier: that is this knob's reduction gate.
    vtg = 3; if (isfield(pa,'clvtgt')), vtg = pa.clvtgt; end
    cp.clvtgt = vtg;    % xpnd_q/fold_p receive cp, NOT pa -- the vent target has
                        % to ride on cp or the vent DOF silently vents nowhere.
    cp.nested = nest;   % fold_p reads this so its pressure pickups stay the
                        % EXACT TRANSPOSE of the xpnd_q injections. If these
                        % two ever disagree the operator loses reciprocity.
    for k=2:(n-1)
        L1_p=aflom(k-1)*pa.chsz(1)./cp.abmom(k); L1_c=aflom(k)*pa.chsz(1)./cp.abmom(k);
        L2_p=aflom(k-1)*pa.chsz(2)./cp.abmom(k); L2_c=aflom(k)*pa.chsz(2)./cp.abmom(k);
        L3_p=aflom(k-1)*pa.chsz(3)./cp.abmom(k); L3_c=aflom(k)*pa.chsz(3)./cp.abmom(k);
        L4_p=aflom(k-1)*pa.chsz(4)./cp.abmom(k); L4_c=aflom(k)*pa.chsz(4)./cp.abmom(k);
        a1(k,1) =-L1_p; a3(k,1) =-L1_c;
        a1(k,6) =-L2_p; a3(k,6) =-L2_c;
        a1(k,11)=-L3_p; a3(k,11)=-L3_c;
        a1(k,16)=-L4_p; a3(k,16)=-L4_c;
        if (nest)
            % NESTED chain ST -d1- CL -d3- SS -d2- SV, with CL VENTING to ST.
            % d1 moves from ST<->SV to ST<->CL so the two new compartments lie
            % BETWEEN the scalae rather than beside them. The vent G represents
            % cortilymph communicating with ST perilymph THROUGH the BM, which
            % the sealed-compartment version omitted; that omission loaded the
            % BM (d1 was ~8x too small at every CL size) and forced the revert.
            % G is scaled by mu3 so pa.clvent is dimensionless; clvent=0 gives
            % the sealed nested chain, and large clvent drives p_CL -> p_ST.
            % Index map (row-1)*4+col with 1=ST 2=SS 3=SV 4=CL.
            G = vent*mu3(k);
            % SEALED nested stamp first, then the vent added symmetrically, so
            % the vent target is a parameter instead of being baked into ST.
            a2(k,1) = L1_p+L1_c+1;            a2(k,4) = -1;            % ST: fluid + d1(ST-CL)
            a2(k,6) = L2_p+L2_c+mu2(k)+mu3(k);a2(k,7) = -mu2(k);       % SS: fluid + d2 + d3
            a2(k,8) = -mu3(k);
            a2(k,10)= -mu2(k);
            a2(k,11)= L3_p+L3_c+mu2(k);                                % SV: fluid + d2 ONLY
            a2(k,13)= -1;                     a2(k,14)= -mu3(k);
            a2(k,16)= L4_p+L4_c+1+mu3(k);                              % CL: fluid + d1 + d3
            % VENT CL <-> chamber vtg. Symmetric off-diagonals keep the operator
            % reciprocal; entries 12 and 15 are untouched by the sealed stamp so
            % they accumulate from zero.
            iv = (vtg-1)*4+vtg; icv = (vtg-1)*4+4; ivc = 12+vtg;
            a2(k,iv)  = a2(k,iv)  + G;   a2(k,16)  = a2(k,16)  + G;
            a2(k,icv) = a2(k,icv) - G;   a2(k,ivc) = a2(k,ivc) - G;
        else
        m3c = cc*mu3(k);   % d3 fluid coupling, scaled by the reduction knob
        a2(k,1) = L1_p+L1_c+1;                a2(k,3) = -1;            % ST: fluid + d1
        a2(k,6) = L2_p+L2_c+mu2(k)+m3c;       a2(k,7) = -mu2(k);       % SS: fluid + d2 + d3
        a2(k,8) = -m3c;
        a2(k,9) = -1;                         a2(k,10)= -mu2(k);
        a2(k,11)= L3_p+L3_c+1+mu2(k);                                  % SV: fluid + d1 + d2
        a2(k,14)= -m3c;                       a2(k,16)= L4_p+L4_c+m3c; % CL: fluid + d3 only
        end
    end
    % Basal boundary: the stapes port acts on ST and SV only (CL and SS have no
    % direct stapes path), mirroring the 3-chamber treatment.
    L1_c=aflom(1)*pa.chsz(1)./cp.abmom(1); L2_c=aflom(1)*pa.chsz(2)./cp.abmom(1);
    L3_c=aflom(1)*pa.chsz(3)./cp.abmom(1); L4_c=aflom(1)*pa.chsz(4)./cp.abmom(1);
    a3(1,1)=-L1_c; a3(1,6)=-L2_c; a3(1,11)=-L3_c; a3(1,16)=-L4_c;
    bcx=1; if (isfield(pa,'bcx')), bcx=pa.bcx; end
    if (nest)
        G1 = vent*mu3(1);
        a2(1,1) = L1_c+1+2*cp.alfx*bcx;     a2(1,4) = -1;
        a2(1,6) = L2_c+mu2(1)+mu3(1);       a2(1,7) = -mu2(1);  a2(1,8) = -mu3(1);
        a2(1,10)= -mu2(1);
        a2(1,11)= L3_c+mu2(1)+2*cp.alfx*bcx;      % stapes still drives SV
        a2(1,13)= -1;                       a2(1,14)= -mu3(1);
        a2(1,16)= L4_c+1+mu3(1);
        iv1 = (vtg-1)*4+vtg; icv1 = (vtg-1)*4+4; ivc1 = 12+vtg;   % same vent, basal row
        a2(1,iv1)  = a2(1,iv1)  + G1;  a2(1,16)   = a2(1,16)   + G1;
        a2(1,icv1) = a2(1,icv1) - G1;  a2(1,ivc1) = a2(1,ivc1) - G1;
    else
    a2(1,1) = L1_c+1+2*cp.alfx*bcx;       a2(1,3) = -1;
    m3c1 = cc*mu3(1);
    a2(1,6) = L2_c+mu2(1)+m3c1;           a2(1,7) = -mu2(1);  a2(1,8) = -m3c1;
    a2(1,9) = -1;                         a2(1,10)= -mu2(1);
    a2(1,11)= L3_c+1+mu2(1)+2*cp.alfx*bcx;
    a2(1,14)= -m3c1;                      a2(1,16)= L4_c+m3c1;
    end
    % Apical boundary: zero pressure gradient in every chamber
    a2(n,:)=0;
    a1(n,1) =-1; a2(n,1) =1;
    a1(n,6) =-1; a2(n,6) =1;
    a1(n,11)=-1; a2(n,11)=1;
    a1(n,16)=-1; a2(n,16)=1;
end
if (m==2), a2(:,(1:m)+(mm-m)) = 1; end
% RAW FLUID BLOCKS, exposed 2026-07-29 so macro_shadow can check them.
% a2's off-diagonal part is a SECOND copy of the chamber-coupling topology that
% macro_couple already holds as Dq/mu/Df -- verified at m=3 that
%     a2(k) == diag(L_p + L_c) + Dq*diag(mu(k,:))*Df
% Two sources of truth for one topology is why pa.rlsplit is ill-posed: it
% changes Df/Dq and this block does not follow. Exposing them is the first step
% to deriving the block instead of writing it out.
cp.a1 = a1; cp.a2 = a2; cp.a3 = a3;
cp.aa=xpnd_a(a1,a2,a3,m,n);
% MACRO/MICRO INTERFACE, built once here and carried on cp because xpnd_q and
% fold_p receive cp, not pa. Everything macro_couple reads (m1/m2/m5/clvm/clvk,
% nested, clcouple, clvtgt) exists by this point: imped() ran at the top of
% cochlea and the chamber stamps are complete.
% (moved earlier -- see the macro_couple call above the fluid operator)
end

% xpnd_a MOVED here from tdm26 (Stage 4): it assembles the block-tridiagonal
% fluid operator from the per-place chamber stamps, which is macromechanics.
% tdm26 no longer referenced it after the extraction -- it would have failed at
% RUNTIME, not parse, because MATLAB local functions are file-scoped.
function aa=xpnd_a(a1,a2,a3,m,n)
% FIXED 'xpnd_a' for m>=1
nm = n*m; nd = 1+2*m;
ad = zeros(nm,nd); dd = zeros(nm,nd);
for j=1:m
    jj = (1:m)+(j-1)*m;
    kk_vec = (j:m:nm)';
    AAA_j = [a1(:, jj), a2(:, jj), a3(:, jj)];
    src_indices = (1:nd) + (j-1); % This was the logic.
    % For m=2, nd=5. j=1 -> 1:5. j=2 -> 2:6.
    % AAA_j has 6 columns. This fits.

    ad(kk_vec, :) = AAA_j(:, src_indices);
end

dd(:, 1+m) = ad(:, 1+m);
for j=1:m
    idx_low = 1+m-j; idx_high = 1+m+j;
    dd(1:(nm-j), idx_low) = ad((1+j):nm, idx_low);
    dd((1+j):nm, idx_high) = ad(1:(nm-j), idx_high);
end
aa = spdiags(dd,-m:m,nm,nm);
end
