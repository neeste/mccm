function R = stiff_diag(nch, boosts, basalfrac)
% STIFF_DIAG  Local partition stiffness / stability diagnostic for tdm26's
%   explicit predictor-corrector integrator (Month-2 scoping).
%
%   For each cochlear place it builds the local 2-DOF partition operator
%   (BM + TM/HB), from the real cp coefficients, as a 4x4 state-space matrix J
%   (state [d1 d2 v1 v2]), and reports:
%     - CONTINUOUS stability: max real part of eig(J) over places (>0 => the
%       model itself is physically unstable, i.e. spontaneous oscillation).
%     - NUMERICAL stability: spectral radius rho(G) of ONE step of the actual
%       predictor-corrector trapezoidal scheme (>1 => the integrator diverges).
%   swept over a basal gampro boost = 1 + boost*ramp (ramp = raised-cosine over
%   the basal fraction, matching abr_regress 'gampro').
%
%   The decisive question: at the observed gampro divergence (nch=1 boost ~0.75)
%   is the continuous system still in the LHP (rho>1 purely NUMERICAL -> an
%   implicit scheme fixes it) or already in the RHP (PHYSICAL -> implicit won't
%   help)?  Fluid mass-loading is omitted, so this LOCAL estimate is a
%   conservative (pessimistic) bound on the boundary.
%
%   Usage:  R = stiff_diag(1, [0 .25 .5 .75 1 1.5], 0.35)

if (nargin<1), nch=1; end
if (nargin<2), boosts=[0 0.25 0.5 0.75 1.0 1.5]; end
if (nargin<3), basalfrac=0.35; end
if (nch>=3), error('stiff_diag: nch=3 (m==3) partition operator not yet implemented'); end

evalc('s = tdm26(0, nch);');       % quiet; get real per-place coefficients
cp = s.cp; pa = s.pa; dt = pa.dt; n = pa.n;
if (ishandle(1)), close(1); end

ramp = 0.5*(1+cos(pi*min(((0:n-1)')/(n-1)/basalfrac, 1)));  % 1 at base -> 0 apically
valid = (cp.m1>0) & (cp.m2>0);
valid(n) = false;                  % skip helicotrema

fprintf('\nstiff_diag nch=%d  dt=%.1e s  n=%d  basalfrac=%.2f\n', nch, dt, n, basalfrac);
% context: fastest local resonance
w1 = sqrt(cp.k1(valid)./cp.m1(valid));            % rad/s
[wmx,~] = max(w1);
fprintf('fastest BM resonance ~ %.1f kHz  (dt*w = %.3f)\n', wmx/2/pi/1000, dt*wmx);
fprintf(['%6s | CONTINUOUS (physical) | NUMERICAL rho(G)     | dt-refine test\n'], '');
fprintf('%6s | %9s %9s | %9s %9s | %10s\n','boost','maxRe/s','verdict','rho(dt)','rho(dt/4)','conclusion');

R.boosts=boosts; R.contRe=nan(size(boosts)); R.rho=nan(size(boosts));
R.rho4=nan(size(boosts)); R.worstk=nan(size(boosts));
for b = 1:numel(boosts)
    gam = cp.gm .* (1 + boosts(b)*ramp);          % cp.gm is the baseline gain (=1)
    maxRe = -inf; maxRho = 0; maxRho4 = 0; worstk = NaN;
    for k = 1:n
        if (~valid(k)), continue; end
        Ba = local_Ba(cp, gam, k);
        maxRe = max(maxRe, max(real(eig([0 0 1 0; 0 0 0 1; Ba]))));
        rk  = max(abs(eig(pc_onestep(Ba, dt))));
        rk4 = max(abs(eig(pc_onestep(Ba, dt/4))));
        if (rk > maxRho), maxRho = rk; worstk = k; end
        maxRho4 = max(maxRho4, rk4);
    end
    R.contRe(b)=maxRe; R.rho(b)=maxRho; R.rho4(b)=maxRho4; R.worstk(b)=worstk;
    % dt-refine test: if rho(dt/4) also >1, shrinking dt does NOT remove it => PHYSICAL
    if (maxRho<=1), concl='stable';
    elseif (maxRho4<=1), concl='NUMERICAL (dt helps)';
    else, concl='PHYSICAL (dt cannot fix)'; end
    fprintf('%6.2f | %9.3g %9s | %9.6f %9.6f | %s\n', boosts(b), maxRe, ...
        tern(maxRe>0,'UNSTABLE','stable'), maxRho, maxRho4, concl);
end
end

% ---- local 2-DOF acceleration operator (m<3), a = Ba*[d1 d2 v1 v2]' ----
function Ba = local_Ba(cp, gam, k)
k1=cp.k1(k);r1=cp.r1(k);m1=cp.m1(k); k2=cp.k2(k);r2=cp.r2(k);m2=cp.m2(k);
k3=cp.k3(k);r3=cp.r3(k); k4=cp.k4(k);r4=cp.r4(k); gh=cp.gh(k); g=gam(k);
Ba = zeros(2,4);
Ba(1,:) = -(1/m1)*[ k1+gh*k3-g*k4,  -gh*k3+g*k4,  r1+gh*r3-g*r4,  -gh*r3+g*r4 ];
Ba(2,:) = -(1/m2)*[ -k3,             k2+k3,        -r3,            r2+r3       ];
end

% ---- one step of the actual predictor-corrector trapezoidal scheme, as G ----
function G = pc_onestep(Ba, dt)
G = zeros(4);
for c = 1:4
    y = zeros(4,1); y(c)=1;
    a  = Ba*y;                              % accel at current
    vp = y(3:4) + a*dt;                     % predictor (forward Euler v)
    dp = y(1:2) + (y(3:4)+vp)*dt/2;
    an = Ba*[dp; vp];                       % accel at predicted state
    vc = y(3:4) + (a+an)*dt/2;              % corrector (trapezoidal)
    dc = y(1:2) + (y(3:4)+vc)*dt/2;
    G(:,c) = [dc; vc];
end
end

function f = bm_cf(cp, k)
if (isnan(k)), f=NaN; return; end
f = sqrt(cp.k1(k)/cp.m1(k))/2/pi;
end
function s=tern(c,a,b), if c, s=a; else, s=b; end, end
