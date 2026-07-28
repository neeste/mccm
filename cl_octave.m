% IS THE CL/d3 RESONANCE THE TIP-TAIL MECHANISM?
%
% SN's hypothesis: the tip-tail distinction arises from a CL resonance, which
% would need to sit roughly 0.5 oct BELOW CF. This sweeps that offset and reads
% the CONTRAST column, which is score26's tip-tail measure and the one column
% every previous sweep in this session ignored.
%
% WHY d3 AND NOT d2. Resonance goes as sqrt(k/m) and each coordinate has a rate
% (ke-me)/2 per cm, with x in cm over xl=3.5:
%   d1  28649 Hz at base, rate -1.2683  ->   338 Hz at apex
%   d2  14646 Hz at base, rate -1.7126  ->  36.5 Hz at apex
% d2's rate is STEEPER than d1's, so its offset below CF drifts from 0.97 oct at
% the base to 3.21 oct at the apex. It cannot hold a fixed relationship to CF,
% so it cannot support a constant tip-tail mechanism along the cochlea. Matching
% d1's rate makes the offset constant by construction, which only d3 is free to
% do: k5/m5 are currently seeded IDENTICAL to k2/m2 and flagged provisional in
% modpar26.m, awaiting exactly this.
%
% CONSTRUCTION. Everything is derived from pa's own k1/m1 so the offset stays
% self-consistent if those are ever refit:
%   k5e = m5e + 2*rate1        (match d1's rate -> constant offset)
%   k5o = m5o * (w1o*2^-off)^2 (place the offset at the base)
%
% CAVEAT this sets the d3 PARTITION resonance relative to the d1 PARTITION
% resonance. The coupled CF differs from sqrt(k1/m1), so a fixed offset here is
% not exactly a fixed offset below realized CF. If contrast peaks off-target,
% that mismatch is the first thing to suspect, not the hypothesis.
%
% OPERATING POINT nested + SV vent at 0.02, og 1: the best gain-per-instability
% row found (ratio 2.04) and stable at maxRe 11.9, so contrast is measured
% without instability confounding it.

A_CRV = [0.95 0.05 0.95 0.05];
OFF = [0 0.25 0.5 0.75 1.0 1.5];

b = modpar26(4);
rate1 = (b.k1e - b.m1e)/2;              % per cm, d1 resonance slope
w1o   = sqrt(b.k1o / b.m1o);            % rad/s at x=0
fprintf('\n  d1 base %.0f Hz, rate %.4f/cm -> apex %.0f Hz (%.2f oct)\n', ...
    w1o/(2*pi), rate1, w1o*exp(rate1*b.xl)/(2*pi), -rate1*b.xl/log(2));
fprintf('  d3 as seeded = d2: base %.0f Hz, rate %.4f/cm\n', ...
    sqrt(b.k5o/b.m5o)/(2*pi), (b.k5e-b.m5e)/2);
fprintf('  CONTRAST is the tip-tail column. Higher = sharper tip over tail.\n\n');
fprintf('  d3 offset      | contrast | amp d1  | maxRe     | range mono | k5o       k5e\n');
fprintf('%s\n', repmat('-',1,88));

R = struct('off',{},'S',{});
% reference row: d3 exactly as seeded from d2 (drifting offset)
for i = 0:numel(OFF)
    pa = modpar26(4); pa.chsz = A_CRV;
    pa.nested = 1; pa.clvtgt = 3; pa.clvent = 0.02; pa.ohcgain = 1;
    if (i == 0)
        lbl = 'as seeded (=d2)';
    else
        off = OFF(i);
        pa.k5e = pa.m5e + 2*rate1;                 % track d1
        pa.k5o = pa.m5o * (w1o * 2^(-off))^2;      % place the offset
        lbl = sprintf('%.2f oct below', off);
    end
    try
        S = score26(pa, 'fast', false);
    catch e
        fprintf('  %-14s | FAILED: %s\n', lbl, e.message); continue
    end
    fprintf('  %-14s | %8.1f | %+6.2f | %+9.1f | %5.2f %-4s | %9.3e %7.4f\n', ...
        lbl, S.contrast, S.amp_gain, S.maxRe, S.bf_range, S.bf_mono, ...
        pa.k5o, pa.k5e);
    R(end+1).off = lbl; R(end).S = S; %#ok<SAGROW>
    save('cl_octave.mat','R');
end
fprintf(['\n  READ THE CONTRAST COLUMN. If it peaks near 0.5 oct, SN''s hypothesis\n' ...
         '  is supported and the offset is a real fitting parameter. If contrast\n' ...
         '  is flat across the whole sweep, the d3 resonance is not what sets\n' ...
         '  tip-tail and the mechanism lies elsewhere. If it rises monotonically\n' ...
         '  to the edge of the sweep, the range needs extending before concluding.\n']);
disp('CL_OCTAVE_DONE');
