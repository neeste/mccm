% RCCHECK -- acceptance tests for the OHC RC pole.
%
% 1. OFF IS BIT-IDENTICAL. No ohctau field must reproduce the committed model
%    exactly, in BOTH solvers. The micro26 rewrite is algebraically identical to
%    the legacy line but NOT bit-for-bit, which is why the off path keeps the
%    original expression verbatim; this checks that it really does.
% 2. tau -> 0 REDUCES TO OFF. A pole far above the audio band must converge back
%    to the instantaneous law. If it does not, the filter is wired wrong.
% 3. THE POLE ACTUALLY DOES SOMETHING, and in the predicted direction: above the
%    corner the active element should become negative DAMPING, which shows up as
%    amplifier gain and stability changing while the CF map does not.
% 4. nimp INDEPENDENCE. The whole reason ohc_rc_step is not inside accel: if the
%    filter were advanced per corrector call, the effective tau would scale with
%    pa.nimp. Changing nimp must NOT change the answer beyond ordinary corrector
%    convergence -- tested against the SAME nimp change with the pole off, so
%    the comparison isolates the filter.
PROJ = '/Users/neely/Library/CloudStorage/OneDrive-FatherFlanagan''sBoysHome/STN/mod/mccm';
addpath(PROJ); cd(PROJ);
L=load('fit_nch3_surface.mat'); pa=L.R.pa;
if (isfield(pa,'hbmode')), pa=rmfield(pa,'hbmode'); end
S=@(p) score26(p,'fast',false);
fprintf('\n=== 1. OFF IS BIT-IDENTICAL ===\n');
R0=fdm26(struct('pa',pa)); s0=S(pa);
fprintf('  fdm26 maperr %.10f | score26 amp %.10f  osc %.10f\n', R0.maperr, s0.amp_gain, s0.maxRe_osc);
fprintf('  (compare against the committed values: maperr 144.02, amp +41.26, osc -180.2)\n');
p0=pa; p0.ohctau=0;                      % explicit zero must also be off
R0b=fdm26(struct('pa',p0)); s0b=S(p0);
d1=[abs(R0b.maperr-R0.maperr) abs(s0b.amp_gain-s0.amp_gain) abs(s0b.maxRe_osc-s0.maxRe_osc)];
fprintf('  ohctau=0 vs absent: maperr %.3e  amp %.3e  osc %.3e  %s\n', d1, ...
        tern(all(d1==0),'IDENTICAL','*** DIFFERS ***'));

fprintf('\n=== 2. tau -> 0 REDUCES TO OFF (fdm26; corner far above the band) ===\n');
for tau=[1e-3 1e-5 1e-7 1e-9]
    p=pa; p.ohctau=tau; R=fdm26(struct('pa',p));
    fprintf('  tau %8.0e (fc %9.3g Hz)  maperr %10.4f   d vs off %10.4f\n', ...
            tau, 1/(2*pi*tau), R.maperr, R.maperr-R0.maperr);
end

fprintf('\n=== 3. DOES THE POLE DO ANYTHING, AND IN WHICH DIRECTION? ===\n');
fprintf('  %10s %10s %10s %10s %10s\n','fc (Hz)','maperr','amp_gain','maxRe_osc','contrast');
taus=[0 1/(2*pi*8000) 1/(2*pi*4000) 1/(2*pi*2000) 1/(2*pi*1000) 1/(2*pi*500)];
for i=1:numel(taus)
    p=pa; if (taus(i)>0), p.ohctau=taus(i); end
    try
        s=S(p);
        fc=Inf; if (taus(i)>0), fc=1/(2*pi*taus(i)); end
        fprintf('  %10.4g %10.2f %+10.2f %10.1f %10.2f\n', ...
                fc, s.maperr, s.amp_gain, s.maxRe_osc, s.contrast);
    catch e
        fprintf('  tau %8.2e FAILED: %s\n', taus(i), e.message(1:min(50,end)));
    end
end

fprintf('\n=== 4. nimp INDEPENDENCE (the reason ohc_rc_step is not inside accel) ===\n');
pon=pa; pon.ohctau=1/(2*pi*1000);
for np=[2 8]
    pa_n=pa;  pa_n.nimp=np;  sn_off=S(pa_n);
    po_n=pon; po_n.nimp=np;  sn_on =S(po_n);
    fprintf('  nimp %2d | OFF amp %+8.3f osc %9.2f | ON amp %+8.3f osc %9.2f\n', ...
            np, sn_off.amp_gain, sn_off.maxRe_osc, sn_on.amp_gain, sn_on.maxRe_osc);
end
fprintf('  READ: the ON rows must move with nimp no more than the OFF rows do.\n');
fprintf('  If ON drifts and OFF does not, the filter is being advanced per corrector.\n');
disp('RCCHECK_DONE');
function s=tern(c,a,b), if c, s=a; else, s=b; end, end
