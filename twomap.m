% DOES THE MODEL PRODUCE A SECOND FREQUENCY MAP?
%
% Allen & Fahey (1993, JASA 94:809-816) report evidence for a second cochlear
% frequency map correlating distortion-product and neural tuning measurements.
% The load-bearing word is MAP. A second RESONANCE only has to exist; a second
% MAP has to be PARALLEL to the first, so that the octave offset is roughly
% constant with place and frequency-to-place is well defined for both.
%
% THE MEASUREMENT. score26's local_tip already computes Am3, the d3 transfer at
% every place, and then never uses it: the BF loop reads only Am (d1). This
% extracts both, place by place, and reports
%   offset(x) = log2( bf1(x) / bf3(x) )
% The MEDIAN offset says where the second map sits. The SPREAD says whether it
% is a map at all. Small spread = parallel = a genuine second map. Large spread
% = a second resonance that drifts = not a map, and the construction fails.
%
% WHY THIS CAN FALSIFY THE CONSTRUCTION. Setting k5e = m5e + 2*rate1 makes the
% d3 PARTITION resonance track d1's by construction. It does NOT follow that the
% REALIZED, fluid-coupled bf3 tracks bf1: coupling can bend either map. The
% as-seeded row (k5 = k2, rate -1.7126 against d1's -1.2683) should show a
% clearly drifting offset. If the rate-matched rows do not show a markedly
% flatter offset, then partition-level rate matching does not survive coupling
% and the parameterization has to be rethought.
%
% METHOD IS COPIED FROM local_tip DELIBERATELY, not reinvented: same f-weighting
% (Am.*f), same 0.15-18 kHz band, same interior-only XLO/XHI restriction. Past
% validation harnesses in this project failed by reducing over frequency instead
% of place, so the validated recipe is reused verbatim.
%
% OPERATING POINT nested, vent 3.0 to SS, ohcgain 1. SN suggested the CL->SS
% vent and it beat the SV vent decisively: amp +55.08 at maxRe +3.9, a ratio of
% 14.1 dB per unit maxRe against m=3b's 4.20 and the SV vent's 3.86, and it
% lands in the 40-60 band at ohcgain 1 with no tuning. My recorded prediction
% that CL->SS would NOT restore the amplifier was wrong: merging CL with SS
% still reaches SV through d2, so ST -d1- (CL+SS) -d2- SV is a complete ST<->SV
% route, just two partitions in series rather than a direct short, and that
% series path is far better damped.

A_CRV = [0.95 0.05 0.95 0.05];
XLO = 0.05; XHI = 0.85;
ISVFRAC = (1391:-10:11)/1401;

b = modpar26(4);
rate1 = (b.k1e - b.m1e)/2;
w1o   = sqrt(b.k1o / b.m1o);

cfg = { 'as seeded (=d2)', NaN
        'rate-matched 0.50 oct', 0.50
        'rate-matched 1.00 oct', 1.00 };

fprintf('\n  d1 partition rate %.4f/cm | d3 as seeded %.4f/cm\n', rate1, (b.k5e-b.m5e)/2);
fprintf('  A SECOND MAP needs the offset spread to be SMALL.\n\n');
fprintf('  config                | median off | spread | IQR  | bf1 rng | bf3 rng | n\n');
fprintf('%s\n', repmat('-',1,80));

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
    if (~isfield(S,'d3') || isempty(S.d3) || all(S.d3(:)==0))
        fprintf('  %-21s | NO d3 SIGNAL\n', cfg{c,1}); continue
    end
    nf = numel(S.f); f = S.f(:); P = fft(S.ped); P = P(1:nf); np = size(S.d1,2);
    xp = pa.isv(:)'/pa.n; inr = xp>=XLO & xp<=XHI;
    bf1 = nan(1,np); bf3 = nan(1,np);
    for i = find(inr)
        D1 = fft(S.d1(:,i)); a1 = abs(D1(1:nf)./max(abs(P),eps)).*f;
        D3 = fft(S.d3(:,i)); a3 = abs(D3(1:nf)./max(abs(P),eps)).*f;
        a1(~(f>0.15&f<18)) = 0; a3(~(f>0.15&f<18)) = 0;
        [q1,i1] = max(a1); [q3,i3] = max(a3);
        if (q1>0), bf1(i) = f(i1); end
        if (q3>0), bf3(i) = f(i3); end
    end
    ok = isfinite(bf1) & isfinite(bf3) & bf1>0 & bf3>0;
    if (sum(ok) < 5), fprintf('  %-21s | too few valid places (%d)\n', cfg{c,1}, sum(ok)); continue; end
    off = log2(bf1(ok)./bf3(ok));
    r1 = log2(max(bf1(ok))/min(bf1(ok)));
    r3 = log2(max(bf3(ok))/min(bf3(ok)));
    fprintf('  %-21s | %+10.3f | %6.3f | %5.3f | %7.2f | %7.2f | %3d\n', ...
        cfg{c,1}, median(off), std(off), iqr(off), r1, r3, sum(ok));
    save(sprintf('twomap_%d.mat',c), 'off', 'bf1', 'bf3', 'pa');
end
fprintf(['\n  READ THE SPREAD. If the rate-matched rows show a much smaller spread\n' ...
         '  than the as-seeded row, partition-level rate matching survives fluid\n' ...
         '  coupling and the model genuinely produces a second MAP, which is the\n' ...
         '  Allen-Fahey observable. The median offset is then a fittable quantity\n' ...
         '  to set against neural and DP data rather than guessed. If all three\n' ...
         '  spreads are comparable, coupling dominates the parameterization and\n' ...
         '  the offset cannot be set through k5 alone.\n']);
disp('TWOMAP_DONE');
