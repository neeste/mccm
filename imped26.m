function cp = imped26(pa)
%IMPED26  Place-dependent partition impedance parameters.
%
% Extracted VERBATIM from tdm26/imped (Stage 4). Belongs with the
% MICROmechanics conceptually -- k1..r5, m1/m2/m5, the vent elements -- but is
% built once during macro setup because the chamber stamps need the mass ratios.

n = pa.n; dx = pa.xl / (n - 1);
x = ((1:n)'-1) * dx; x = x.*(1+(pa.xtap*x).^pa.xtex); q = x .^ 2;
cp.k1 = pa.k1o * exp(pa.k1e * x + pa.k1q * q);
if isfield(pa, 'rough_amp')                       % coherent-reflection roughness (ported from fdm26/imped)
    sd=42; if (isfield(pa,'rough_seed')), sd=pa.rough_seed; end
    rng(sd);                                      % reproducible place-fixed roughness (seed selects the realization)
    rr = rand(size(cp.k1));                       % realization (full length -> seed-consistent)
    if (isfield(pa,'rough_xlo'))                  % optional basal cutoff: roughness only apical of frac*L
        rr(( 1:numel(cp.k1))' < round(pa.rough_xlo*numel(cp.k1))) = 0;
    end
    if (isfield(pa,'rough_fc') && pa.rough_fc>0)  % optional per-CF band (DIAGNOSTIC probe): roughness
        fcf = sqrt((pa.k1o*exp(pa.k1e*x+pa.k1q*q)) ...   % only within +-rough_foct octaves of the
                 ./(pa.m1o*exp(pa.m1e*x+pa.m1q*q)))/(2*pi*1000);  % stimulus CF place (BM-resonance map, kHz)
        woct=0.75; if (isfield(pa,'rough_foct')), woct=pa.rough_foct; end
        rr(abs(log2(fcf/pa.rough_fc)) > woct) = 0;
    end
    cp.k1 = cp.k1 .* (1 + pa.rough_amp * rr);     % opt-in: absent pa.rough_amp -> smooth (unchanged)
end
cp.r1 = pa.r1o * exp(pa.r1e * x + pa.r1q * q);
cp.m1 = pa.m1o * exp(pa.m1e * x + pa.m1q * q);
cp.k2 = pa.k2o * exp(pa.k2e * x + pa.k2q * q);
cp.r2 = pa.r2o * exp(pa.r2e * x + pa.r2q * q);
cp.m2 = pa.m2o * exp(pa.m2e * x + pa.m2q * q);
cp.k3 = pa.k3o * exp(pa.k3e * x + pa.k3q * q);
cp.r3 = pa.r3o * exp(pa.r3e * x + pa.r3q * q);
cp.k4 = pa.k4o * exp(pa.k4e * x + pa.k4q * q);
cp.r4 = pa.r4o * exp(pa.r4e * x + pa.r4q * q);
if (isfield(pa,'k5o'))   % DOF-3 (OC height / cortilymph pump) -- 4-chamber only
    cp.k5 = pa.k5o * exp(pa.k5e * x + pa.k5q * q);
    cp.r5 = pa.r5o * exp(pa.r5e * x + pa.r5q * q);
    cp.m5 = pa.m5o * exp(pa.m5e * x + pa.m5q * q);
end
% ---- VENT AS A RESONANT ELEMENT (DOF 4) ---------------------------------
% The vent was ALREADY inertial: its a2 stamp G = clvent*mu3 = clvent*m1/m5 is
% exactly the coefficient m1/mv of a DOF whose mass is mv = m5/clvent, algebra-
% ically eliminated because it had no stiffness or damping. That is why the CL
% resonance could not be placed: with a pure inertance, ONE number sets both the
% coupling strength and the resonance, so tightening the coupling necessarily
% moves the resonance. twomap2 measured the consequence -- the realized d1->d3
% transfer peak sat at CF regardless of k5 (medians -0.060 vs -0.114 oct for
% partition resonances a full octave apart, indistinguishable at IQR ~0.5).
%
% Giving the vent a stiffness adds the second element, so the resonance
% sqrt(clvk/mv) becomes placeable INDEPENDENTLY of the coupling strength. This
% is a Helmholtz resonator on the cortilymph space: mv is the channel inertance,
% clvk the compliance of the space, clvr the channel resistance.
%
% REDUCTION GATE clvk = clvr = 0 leaves the vent a pure inertance and must
% reproduce the previous behaviour EXACTLY. The a2 stamp is unchanged in that
% limit (still G), and s4 is identically zero, so nothing is injected.
if (isfield(pa,'clvent') && pa.clvent > 0 && isfield(cp,'m5'))
    vv = pa.clvent;
    cp.clvm = cp.m5 ./ vv;                       % channel inertance, mv = m5/clvent
    kv = 0; if (isfield(pa,'clvk')), kv = pa.clvk; end
    rv = 0; if (isfield(pa,'clvr')), rv = pa.clvr; end
    % pa.clvoct places the CL resonance a fixed number of octaves BELOW the
    % local BM resonance sqrt(k1/m1), and takes precedence over a raw clvk. It
    % is resolved here, not in modpar26, because the required stiffness is
    % place-dependent (it tracks k1/m1/m5) so a scalar in the parameter file
    % would be right at exactly one place. fdm26/imped resolves it identically.
    if (isfield(pa,'clvoct') && isfinite(pa.clvoct))
        kv = cp.clvm .* (cp.k1 ./ cp.m1) * 4^(-pa.clvoct);
    end
    % clvk/clvr may be scalars or n-vectors; scalars are broadcast along place.
    cp.clvk = kv .* ones(size(cp.m5));
    cp.clvr = rv .* ones(size(cp.m5));
end
cp.gh = pa.gpo * exp(pa.gpe * x + pa.gpq * q);
cp.ac = pa.aco * exp(pa.ace * x + pa.acq * q);
cp.bw = pa.bwo * exp(pa.bwe * x + pa.bwq * q);
if (isfield(pa,'z0unif') && pa.z0unif)
    % Uniform-z0 area compensation (ported from fdm26/imped). Retaper scala area
    % ac and BM width bw by the SAME factor exp(sz*x): bw/ac -- hence the coupling
    % aflom./abmom, kappa=sqrt(Zs*Yv), and the CF map / tuning / latency -- is
    % preserved to leading order, while z0=sqrt(Zs/Yv) is flattened in the
    % stiffness tail, killing the spurious frequency-independent BASAL reflection
    % that otherwise dominates the emission. pa.z0slope = MEASURED acoustic-mode z0
    % log-slope (correct for the 3-chamber's 2-DOF coupling); default forces
    % ac*bw ~ k1 (exact only when partition admittance Yb ~ s/k1).
    if (isfield(pa,'z0slope'))
        sz=pa.z0slope; szq=0; if (isfield(pa,'z0slopeq')), szq=pa.z0slopeq; end
        ace=pa.ace+sz; bwe=pa.bwe+sz; acq=pa.acq+szq; bwq=pa.bwq+szq;
    else
        d =pa.bwe-pa.ace; ace=(pa.k1e-d)/2; bwe=(pa.k1e+d)/2;
        dq=pa.bwq-pa.acq; acq=(pa.k1q-dq)/2; bwq=(pa.k1q+dq)/2;
    end
    cp.ac=pa.aco*exp(ace*x+acq*q);
    cp.bw=pa.bwo*exp(bwe*x+bwq*q);
end
cp.gm = pa.gam * pa.gampro;
ihe=n;
cp.k1(ihe) = pa.khe; cp.r1(ihe) = pa.rhe; cp.m1(ihe) = pa.mhe;
cp.k2(ihe) = 0; cp.r2(ihe) = 0; cp.m2(ihe) = 0; cp.k3(ihe) = 0;
cp.r3(ihe) = 0; cp.k4(ihe) = 0; cp.r4(ihe) = 0;
if (isfield(cp,'k5')), cp.k5(ihe) = 0; cp.r5(ihe) = 0; end
end

