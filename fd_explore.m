function fd_explore
% Month-2 group-delay exploration: which dispersion parameter moves the
% frequency-domain forward-latency slope d toward the ABR target 0.413, and at
% what near-criticality cost (tipgain)? One-at-a-time sweeps around baseline.
% d=0.690 too steep; LOWER d = more basal delay = toward target.
base = modpar26(1);
b = fd_grpdelay(base);
fprintf('BASELINE: d=%.3f  tipgain=%.1f dB   (target d=0.413)\n\n', b.d, b.tipgain);
sweeps = { ...
  'k1e',  [-2.40 -2.55 -2.6655 -2.80 -3.00];   % BM stiffness gradient
  'm1e',  [-0.05 -0.10 -0.1884 -0.30 -0.50];   % BM mass gradient
  'r1e',  [-1.00 -1.20 -1.3937 -1.60 -1.90];   % BM damping gradient
  'xtap', [ 0.10  0.18  0.2339  0.30  0.40];    % place-warp strength
  'xtex', [ 3     4     6       8     10  ]};   % place-warp sharpness
flst = logspace(log10(200), log10(20000), 500);
for s = 1:size(sweeps,1)
    fld = sweeps{s,1}; vals = sweeps{s,2};
    fprintf('== %-5s (base %.4g) ==   %8s %8s\n', fld, base.(fld), 'slope_d','tipgain');
    for v = vals
        pa = base; pa.(fld) = v;
        R = fd_grpdelay(pa, flst);
        mark = ''; if (abs(R.d-0.413) < abs(b.d-0.413)-0.02), mark = ' <-toward target'; end
        fprintf('   %-5s = %+7.4g       %8.3f %8.1f%s\n', fld, v, R.d, R.tipgain, mark);
    end
    fprintf('\n');
end
end
