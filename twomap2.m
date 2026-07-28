% THE SECOND MAP, MEASURED AS THE d1 -> d3 TRANSFER.
%
% WHY twomap WAS THE WRONG MEASUREMENT. It took the GLOBAL argmax of |d1| and of
% |d3| separately and compared them. Result: median offset exactly 0.000 with
% IQR 0.000 at 112 places, in all three configurations. That is expected almost
% by construction -- in a travelling-wave cochlea every partition DOF at a place
% is driven by the same wave and peaks at the local CF, so the global maxima
% coincide. A second map present as a SUBDOMINANT feature (a shoulder or
% secondary peak on d3) is invisible to a global argmax. The negative result
% says the test was blunt, not that the hypothesis is wrong.
%
% THE RIGHT MEASUREMENT. Divide out the common envelope:
%     H32(f) = |d3(f)| / |d1(f)|
% This is the transfer from BM to OC height, and its peak is whatever d3 adds on
% top of the BM. That is much closer to what Allen & Fahey (1993) expose, since
% DP and neural tuning are read RELATIVE to the BM response rather than in
% absolute terms.
%
% GUARD, and it matters. H32 blows up wherever d1 -> 0, which happens far from
% CF and would put the peak at a numerical artifact rather than a resonance.
% Every frequency where |d1| falls below DFLOOR of its own peak at that place is
% excluded. Without this the test finds the edge of the passband every time.
%
% READ: offset = log2(bf1 / bf_H32) in octaves, positive meaning the transfer
% peaks BELOW CF. SN's hypothesis puts it near 0.5. The IQR is the map test --
% a genuine second MAP needs the offset roughly constant across place, not just
% nonzero somewhere.

A_CRV = [0.95 0.05 0.95 0.05];
XLO = 0.05; XHI = 0.85; DFLOOR = 0.05;
ISVFRAC = (1391:-10:11)/1401;

b = modpar26(4);
rate1 = (b.k1e - b.m1e)/2;
w1o   = sqrt(b.k1o / b.m1o);

cfg = { 'as seeded (=d2)', NaN
        'rate-matched 0.50 oct', 0.50
        'rate-matched 1.00 oct', 1.00 };

fprintf('\n  H32 = |d3|/|d1|, peak measured only where |d1| >= %.2f of its peak.\n', DFLOOR);
fprintf('  offset = log2(bf1/bf_H32), POSITIVE = transfer peaks BELOW CF.\n');
fprintf('  SN hypothesis: near +0.5 oct, and CONSTANT across place (small IQR).\n\n');
fprintf('  config                | median off | IQR   | spread | frac>0 | n\n');
fprintf('%s\n', repmat('-',1,72));

for c = 1:size(cfg,1)
    pa = modpar26(4); pa.chsz = A_CRV;
    pa.nested = 1; pa.clvtgt = 2; pa.clvent = 3.0; pa.ohcgain = 1;
    if (isfinite(cfg{c,2}))
        pa.k5e = pa.m5e + 2*rate1;
        pa.k5o = pa.m5o * (w1o * 2^(-cfg{c,2}))^2;
    end
    pa.isv = unique(max(1, min(pa.n, round(ISVFRAC*pa.n))), 'stable');
    try
        evalc('S = tdm26(0,pa,0,0);');
    catch e
        fprintf('  %-21s | THREW: %s\n', cfg{c,1}, e.message); continue
    end
    if (any(~isfinite(S.d1(:)))), fprintf('  %-21s | DIVERGED\n', cfg{c,1}); continue; end
    nf = numel(S.f); f = S.f(:); P = fft(S.ped); P = P(1:nf); np = size(S.d1,2);
    xp = pa.isv(:)'/pa.n; inr = xp>=XLO & xp<=XHI;
    off = nan(1,np);
    for i = find(inr)
        A1 = abs(fft(S.d1(:,i))); A1 = A1(1:nf);
        A3 = abs(fft(S.d3(:,i))); A3 = A3(1:nf);
        bnd = f>0.15 & f<18;
        a1w = A1.*f; a1w(~bnd) = 0;
        [p1,i1] = max(a1w); if (p1<=0), continue; end
        ok = bnd & (A1 >= DFLOOR*max(A1(bnd)));    % the guard
        if (sum(ok) < 5), continue; end
        H = nan(nf,1); H(ok) = A3(ok)./A1(ok);
        [pH,iH] = max(H); if (~isfinite(pH)), continue; end
        off(i) = log2(f(i1)/f(iH));
    end
    v = off(isfinite(off));
    if (numel(v) < 5), fprintf('  %-21s | too few valid places (%d)\n', cfg{c,1}, numel(v)); continue; end
    fprintf('  %-21s | %+10.3f | %5.3f | %6.3f | %5.2f | %3d\n', ...
        cfg{c,1}, median(v), iqr(v), std(v), mean(v>0), numel(v));
    save(sprintf('twomap2_%d.mat',c), 'off', 'pa');
end
fprintf(['\n  A SECOND MAP needs median offset clearly nonzero AND a small IQR.\n' ...
         '  Nonzero median with large IQR = a second resonance that drifts, which\n' ...
         '  is not a map. Median near zero = d3 adds no resonance of its own on\n' ...
         '  top of the BM, and the tip-tail mechanism is not the CL resonance.\n']);
disp('TWOMAP2_DONE');
