% GPTEST -- round-trip and no-op checks for the gampro-grade parameter vector.
%
% The no-op check is the one that matters: grade 0 must reproduce ones(n,1)
% EXACTLY (not approximately), so every fit saved before the parameter existed
% still means what it meant, and so an old-length pv still loads. Exercises the
% shipped mapping via parfit26('handles') rather than a re-implementation.
H = parfit26('handles');
fprintf('\n== gampro-grade pv checks ==\n');
allok = true;
for nch = [1 2 3]
    pa  = modpar26(nch);
    pv  = H.getpar(pa);
    nc  = numel(pa.chsz); gpi = 30 + nc + 1;   % gampro grade; hbsc is gpi+1
    lenok = (numel(pv) == gpi + 1);

    pa2 = H.setpar(pa, pv, []);          % round trip at the default grade (0)
    rt2 = H.getpar(pa2);
    d30 = max(abs(pv(1:30) - rt2(1:30)));
    dch = max(abs(pa.chsz - pa2.chsz));
    dgp = max(abs(pa2.gampro - 1));      % must be EXACTLY 0

    pa3  = H.setpar(pa, pv(1:30+nc), []); % pre-2026-08-01 pv length
    dgp3 = max(abs(pa3.gampro - 1));

    pvg = pv; pvg(gpi) = -0.20;           % a real grade must bite, and round trip
    pa4 = H.setpar(pa, pvg, []);
    rt4 = H.getpar(pa4);
    rng4 = [min(pa4.gampro) max(pa4.gampro)];
    dgr  = abs(rt4(gpi) + 0.20);

    % hbsc: default preserved, round trips, and clamped at 0 (a negative slope
    % would make gain GROW with level -- expansive, not compressive)
    pvh = pv; pvh(gpi+1) = 0.64;  pa5 = H.setpar(pa, pvh, []);
    pvn = pv; pvn(gpi+1) = -1.0;  pa6 = H.setpar(pa, pvn, []);
    hbok = abs(pv(gpi+1)-0.04) < 1e-15 && abs(pa5.hbsc-0.64) < 1e-15 && pa6.hbsc == 0;

    ok = lenok && d30 < 1e-15 && dch < 1e-15 && dgp == 0 && dgp3 == 0 ...
         && dgr < 1e-15 && rng4(1) < 1 && rng4(2) > 1 && hbok;
    allok = allok && ok;
    fprintf('  nch=%d nc=%d gpi=%2d | d30 %.1e dchsz %.1e | grade0 max|gp-1| %.1e | oldpv %.1e', ...
            nch, nc, gpi, d30, dch, dgp, dgp3);
    fprintf(' | g=-0.20 range [%.4f %.4f] rt %.1e | %s\n', ...
            rng4(1), rng4(2), dgr, tern(ok,'PASS','FAIL'));
end
fprintf('\n  %s\n', tern(allok, 'ALL PASS', 'FAILURES ABOVE'));
disp('GPTEST_DONE');
function s = tern(c,a,b), if c, s=a; else, s=b; end, end
