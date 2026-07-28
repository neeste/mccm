% Where does the raw-surface misfit live?  Compare the model's 16 latencies
% (+Delta) against the 1988/2013 target band, cell by cell.  No model run.
L=load('tmp_surface_check.mat'); R=L.R; m=R.mf;
f=m.f(:); slv=m.slv(:);
[F,I]=ndgrid(f, slv/100);
t88 = 12.90*(5.00.^(-I)).*(F.^(-0.413));
t13 = 12.63*(5.34.^(-I)).*(F.^(-0.390));
tlo=min(t88,t13); thi=max(t88,t13);
D=R.delta; T=m.lat+D;
V=max(0,tlo-T)+max(0,T-thi);
fprintf('Delta=%.2f ms   band-RMS=%.3f ms\n\n', D, sqrt(mean(V(:).^2)));
fprintf(' f(kHz) lvl   model  +Delta   band_lo band_hi   viol\n');
for a=1:numel(f)
  for b=1:numel(slv)
    fprintf('%6.2f %3.0f  %6.2f  %6.2f   %6.2f  %6.2f  %6.2f\n', ...
        f(a), slv(b), m.lat(a,b), T(a,b), tlo(a,b), thi(a,b), V(a,b));
  end
end
fprintf('\nRMS violation by frequency:\n');
for a=1:numel(f), fprintf('  %.2f kHz : %.3f ms\n', f(a), sqrt(mean(V(a,:).^2))); end
fprintf('RMS violation by level:\n');
for b=1:numel(slv), fprintf('  %3.0f dB  : %.3f ms\n', slv(b), sqrt(mean(V(:,b).^2))); end
disp('SURFACE_DIAG_DONE');
