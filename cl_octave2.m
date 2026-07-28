% TIP-TAIL vs d3 OFFSET, WITH THE AMPLIFIER HELD FIXED.
%
% WHY cl_octave DID NOT TEST THE HYPOTHESIS. Placing the d3 resonance means
% changing k5o, and k5o also sets loop gain. That sweep moved BOTH: amp ran from
% +24 to +238 dB and maxRe from +11.9 to +291.8. Contrast measured at maxRe
% 291.8 is not comparable to contrast at 11.9, so the contrast ordering it
% produced (minimum at 0.5 oct, maximum at 0.0) may be entirely a gain artifact.
% Same confounding that invalidated the first nested test.
%
% THE CONTROL. At each offset, bisect ohcgain until amp d1 hits AMPTGT. Then the
% only thing differing between rows is where the d3 resonance sits. maxRe is
% reported too: if it also equalizes the comparison is clean; if not, proximity
% to criticality still differs and the control is only partial.
%
% OPERATING POINT UPDATED to the CL->SS vent, which vent_ss showed is the better
% construction: amp +52 to +55 dB at ohcgain 1 with maxRe 3.9-7.7, a ratio up to
% 14.1 dB per unit maxRe against m=3b's 4.20. The earlier SV vent needed
% ohcgain 0.85 to reach +44.67 at maxRe 11.6 (ratio 3.86).
%
% BUG FIXED local functions declared in a SCRIPT do not share the script's
% workspace, unlike nested functions inside a function file. The previous
% version put the pa construction in a local mk() and it could not see A_CRV,
% rate1 or w1o. Inlined here.
%
% HYPOTHESIS UNDER TEST (SN): the CL resonance is the tip-tail mechanism, which
% would put it near 0.5 oct below CF. Supported if contrast peaks near 0.5 oct
% once gain is controlled. Refuted if contrast is flat or peaks far from there.

A_CRV = [0.95 0.05 0.95 0.05];
AMPTGT = 45; NIT = 6; OGLO = 0.02; OGHI = 1.6;
OFF = [0 0.25 0.5 0.75 1.0 1.5];

b = modpar26(4);
rate1 = (b.k1e - b.m1e)/2;
w1o   = sqrt(b.k1o / b.m1o);

fprintf('\n  AMPLIFIER HELD AT %g dB by bisecting ohcgain. Only the d3 offset varies.\n', AMPTGT);
fprintf('  Operating point: nested, vent 3.0 -> SS (best gain per unit maxRe).\n');
fprintf('  Uncontrolled sweep said: contrast 32.8 at 0.0 oct, 12.2 at 0.5 oct.\n\n');
fprintf('  d3 offset       | ohcgain | amp d1  | contrast | maxRe    | range mono\n');
fprintf('%s\n', repmat('-',1,78));

% Accumulate and rewrite the whole struct each row. The previous version used
% save(...,'-append'), which REQUIRES the file to already exist and threw after
% row 1 because the launch command had just deleted it.
R = struct('lbl',{},'og',{},'S',{});
ALL = [NaN OFF];
for c = 1:numel(ALL)
    off = ALL(c);
    if (isfinite(off)), lbl = sprintf('%.2f oct below', off); else, lbl = 'as seeded (=d2)'; end
    lo = OGLO; hi = OGHI; Sb = []; ogb = NaN;
    for it = 1:NIT
        og = 0.5*(lo+hi);
        pa = modpar26(4); pa.chsz = A_CRV;
        pa.nested = 1; pa.clvtgt = 2; pa.clvent = 3.0;   % CL -> SS
        if (isfinite(off))
            pa.k5e = pa.m5e + 2*rate1;
            pa.k5o = pa.m5o * (w1o * 2^(-off))^2;
        end
        pa.ohcgain = og;
        try
            S = score26(pa, 'fast', false);
        catch
            hi = og; continue
        end
        if (~isfinite(S.amp_gain) || S.maxRe > 5e3)
            hi = og; continue
        end
        if (isempty(Sb) || abs(S.amp_gain-AMPTGT) < abs(Sb.amp_gain-AMPTGT))
            Sb = S; ogb = og;
        end
        if (S.amp_gain < AMPTGT), lo = og; else, hi = og; end
    end
    if (isempty(Sb))
        fprintf('  %-15s | no usable ohcgain found in [%.2f %.2f]\n', lbl, OGLO, OGHI);
        continue
    end
    fprintf('  %-15s | %7.4f | %+6.2f | %8.1f | %+8.1f | %5.2f %-4s\n', ...
        lbl, ogb, Sb.amp_gain, Sb.contrast, Sb.maxRe, Sb.bf_range, Sb.bf_mono);
    R(end+1).lbl = lbl; R(end).og = ogb; R(end).S = Sb; %#ok<SAGROW>
    save('cl_octave2.mat','R');
end
fprintf(['\n  READ: with amp pinned, does contrast peak near 0.5 oct? Check the amp\n' ...
         '  column first -- rows that missed the target are not comparable and the\n' ...
         '  bisection range may need widening. Check maxRe second: if it varies a\n' ...
         '  lot at fixed amp, proximity to criticality still differs between rows\n' ...
         '  and the control is only partial.\n']);
disp('CL_OCTAVE2_DONE');
