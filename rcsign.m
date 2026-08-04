% RCSIGN -- does the pole actually produce NEGATIVE DAMPING, and what does it do
% to stability? rccheck could not answer either: it read maxRe_osc from
% score26('fast',false), which does not compute it (NaN in every row), and its
% pass/fail line used all(d==0) on a vector containing NaN, so it printed
% "DIFFERS" for two quantities that were both exactly 0.000e+00.
%
% THE SIGN QUESTION IS THE ONE THAT MATTERS. The claim is that Za = Ka/(i*omega)
% (90 deg, negative stiffness) plus one RC pole (another 90 deg above corner)
% gives 180 deg = negative DAMPING. Negative damping means the active term
% contributes a NEGATIVE REAL PART to the partition impedance. That is checkable
% directly from imped, with no fit and no ambiguity -- and it must be checked,
% because rccheck showed amp_gain FALLING (+41.3 -> +25.3), which is what an
% INVERTED sign would also look like.
PROJ = '/Users/neely/Library/CloudStorage/OneDrive-FatherFlanagan''sBoysHome/STN/mod/mccm';
addpath(PROJ); cd(PROJ);
L=load('fit_nch3_surface.mat'); pa=L.R.pa;
if (isfield(pa,'hbmode')), pa=rmfield(pa,'hbmode'); end
tau=1/(2*pi*1000);
x=0.5*pa.xl;                          % mid-cochlea place
fprintf('\n  effective ACTIVE impedance  -gam*za  at x=L/2, tau=%.3e (fc 1000 Hz)\n', tau);
fprintf('  negative Re = energy INJECTED (negative damping); positive Re = dissipation\n\n');
fprintf('  %8s | %14s %14s | %14s %14s\n','f (Hz)','Re off','Im off','Re on','Im on');
for f=[250 500 1000 2000 4000 8000]
    s=2i*pi*f;
    k4=pa.k4o*exp(pa.k4e*x+pa.k4q*x^2); r4=pa.r4o*exp(pa.r4e*x+pa.r4q*x^2);
    za_off=k4/s+r4;
    za_on =za_off/(1+s*tau);
    g=pa.gam;
    aoff=-g*za_off; aon=-g*za_on;
    fprintf('  %8.0f | %14.4g %14.4g | %14.4g %14.4g\n', f, real(aoff), imag(aoff), real(aon), imag(aon));
end
fprintf('\n=== STABILITY with the pole on (coupeig, which score26 fast does not run) ===\n');
for fc=[Inf 4000 1000]
    p=pa; if (isfinite(fc)), p.ohctau=1/(2*pi*fc); end
    try
        evalc('C=tdm26(''coupeig'',struct(''pa'',p));');
        R=fdm26(struct('pa',p));
        fprintf('  fc %8.4g Hz | maxRe %+9.2f | maxRe_osc %+9.2f | maperr %8.2f\n', ...
                fc, C.maxRe, C.maxRe_osc, R.maperr);
    catch e
        fprintf('  fc %8.4g Hz FAILED: %s\n', fc, e.message(1:min(60,end)));
    end
end
disp('RCSIGN_DONE');
