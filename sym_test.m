% IS m=1 REALLY m=2 UNDER A SYMMETRY CONSTRAINT? (SN, 2026-07-25)
% If so the two ARE genuinely nested (unlike m=2 vs m=3, which use different
% force laws), and the m=1 champion's maperr 74.6 must be reachable by m=2 --
% making the m=2 champion (parfit26_lowdamp, 191.3) simply under-optimized.
L=load('refit_c1broad.mat'); pc=L.R.pa;
fprintf('m=1 champion: pa.m=%g  chsz=[%s]  (par_CEL16 default is [1 1] = symmetric)\n', ...
        pc.m, num2str(pc.chsz,'%.4f '));
cfg={ {'m=1 native',              setfield(modpar26(1),'m',1)}, ...
      {'m=2 native',              modpar26(2)}, ...
      {'m=1 CHAMPION',            pc}, ...
      {'m=1 champion params AS m=2', setfield(pc,'m',2)} };
fprintf('\n  config                        | maperr | range mono | contr | amp    | maxRe\n');
fprintf('%s\n',repmat('-',1,84));
for c=1:numel(cfg)
    nm=cfg{c}{1}; pa=cfg{c}{2};
    try
        S=score26(pa,'fast',false);
        fprintf('  %-29s | %6.1f | %5.2f %-4s | %5.1f | %+6.2f | %+7.1f\n', ...
            nm,S.maperr,S.bf_range,S.bf_mono,S.contrast,S.amp_gain,S.maxRe);
    catch e
        fprintf('  %-29s | FAILED: %s\n', nm, e.message);
    end
end
fprintf(['\n  IF m=1 == m=2-with-symmetry, rows 1 and 2 should MATCH (same par_CEL16\n' ...
         '  params, chsz=[1 1]) and row 4 should reproduce row 3 (74.6). Row 4 beating\n' ...
         '  the m=2 champion 191.3 would prove that champion is under-optimized.\n']);
disp('SYM_TEST_DONE');
