% Validate abr_surface_obj on the saved fits (zero compute -- uses stored mf).
% Question: does fitting the FULL 16-point surface give a different picture than
% the summary slope, and is the model's latency surface power-law SHAPED at all?
files={'parfit26_shoulder.mat','parfit26_recip.mat'};
for k=1:numel(files)
    if (~exist(files{k},'file')), fprintf('(%s missing)\n',files{k}); continue; end
    L=load(files{k}); m=L.R.mf;
    [J,D]=abr_surface_obj(m);
    fprintf('\n=== %s ===\n',files{k});
    fprintf('  summary-stat objective saw : slope=%.3f  level_c=%.3f (%.3f %%/dB)\n', ...
            m.slope, m.level_c, 100*(m.level_c^(1/100)-1));
    fprintf('  FULL-surface fit  : b=%.2f  c=%.3f (%.3f %%/dB)  d=%.3f   [n=%d/16 cells]\n', ...
            D.b, D.c, 100*(D.c^(1/100)-1), D.d, D.n);
    fprintf('  targets           : b 11.99-13.27   c 5.00-5.34    d 0.39-0.41\n');
    fprintf('  SHAPE error (rms log-resid) = %.4f   (how power-law-like the surface is)\n', D.resid);
    fprintf('  hinge losses      : d=%.4f  c=%.4f  b=%.4f   -> J_surface=%.4f\n', ...
            D.hd,D.hc,D.hb,J);
    fprintf('  cells excluded as detector failures: %d\n', 16-D.n);
    % where does the surface deviate most?
    R=nan(size(m.lat)); R(D.mask)=D.r;
    fprintf('  per-cell log-residual (rows=f %s kHz, cols=lv %s dB):\n', ...
            mat2str(m.f(:)'), mat2str(m.slv(:)'));
    for i=1:size(R,1)
        fprintf('    ');
        for j=1:size(R,2)
            if (isnan(R(i,j))), fprintf('    xx '); else, fprintf('%+7.3f ',R(i,j)); end
        end
        fprintf('\n');
    end
end
disp('SURFACE_DEMO_DONE');
