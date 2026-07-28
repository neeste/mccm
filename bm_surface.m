% The corrected latency surface: BM-peak measurement (hbmode='bm', verified to
% leave the mechanics bit-identical) on the current best parameters, scored
% against the 1988/2013 target band with a free neural-delay offset Delta.
L=load('parfit26_recip.mat'); pa=L.R.pa; pa.hbmode='bm';
m=abr_metric(pa,false);
if (~m.ok), fprintf('abr_metric failed: %s\n', m.msg); return; end
f=m.f(:); slv=m.slv(:);
[F,I]=ndgrid(f, slv/100);
t88=12.90*(5.00.^(-I)).*(F.^(-0.413));
t13=12.63*(5.34.^(-I)).*(F.^(-0.390));
tlo=min(t88,t13); thi=max(t88,t13);
Lm=m.lat; ok=isfinite(Lm)&Lm>0;
ds=linspace(0,1.5,301); best=inf; D=0;
for k=1:numel(ds)
    V=max(0,tlo-(Lm+ds(k)))+max(0,(Lm+ds(k))-thi);
    r=sqrt(mean(V(ok).^2)); if (r<best), best=r; D=ds(k); end
end
V=max(0,tlo-(Lm+D))+max(0,(Lm+D)-thi);
fprintf('BM-peak surface: Delta=%.2f ms (a=%.2f)  band-RMS=%.3f ms  %d/%d cells\n\n',D,5+D,best,sum(ok(:)),numel(Lm));
fprintf(' f(kHz) lvl   lat_bm  +Delta   band_lo band_hi   viol\n');
for a=1:numel(f)
  for b=1:numel(slv)
    fprintf('%6.2f %3.0f  %6.2f  %6.2f   %6.2f  %6.2f  %6.2f\n',f(a),slv(b),Lm(a,b),Lm(a,b)+D,tlo(a,b),thi(a,b),V(a,b));
  end
end
fprintf('\nper-20dB %% change in latency (model vs data ~28%%/20dB):\n');
for a=1:numel(f)
  fprintf('  %.2f kHz : ',f(a));
  for b=1:numel(slv)-1, fprintf('%6.1f%%',100*(Lm(a,b+1)-Lm(a,b))/Lm(a,b)); end
  fprintf('\n');
end
fprintf('\ndiagnostic power-law: slope=%.3f  level_c=%.2f (%.2f %%/dB)\n', m.slope, m.level_c, 100*(m.level_c^0.01-1));
save('bm_surface.mat','m','D','best');
disp('BM_SURFACE_DONE');
