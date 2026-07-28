% BISECT THE m=4 AMPLIFIER FROM m=3.
%
% m=3 and m=4 share their impedance profile to within ~15%, yet m=3 amplifies
% +81 dB (d1) / +84 dB (d2) and m=4 manages +2.40 / +0.60 / -5.32. The structural
% difference is the third DOF and its couplings. This locates the fault.
%
% HYPOTHESIS, from reading the two force laws side by side:
%   m=3 : s1 = -(k1*d1 + r1*v1 + k_act*d2 + r_act*v2)   BM receives -act
%   m=4 : s1 = -(k1*d1 + r1*v1 - act)                   BM receives +act
% THE BM FORCE SIGN IS INVERTED between them. That explains both earlier
% results: ohcsgn=+1 is dissipative because it is backwards relative to the
% model that achieves +81 dB, and ohcsgn=-1 restores the m=3 sense, which is why
% it injects (ohcP > 0). It then diverges because m=4 ALSO puts +act on d3, a
% second energy path m=3 does not have.
%
% PREDICTION: with ohcsgn=-1 (m=3 BM sense) AND d3 progressively frozen by a
% heavy m5, the d3 branch stops absorbing and destabilizing, and m=4 should
% approach m=3's amplifier while staying stable. If it does, the fault is the
% d3 reaction path and it is localized. If amp stays low at every m5, the sign
% is not the whole story and the fault lies elsewhere.
%
% chsz(4) stays at its native 0.05: setting it to 0 algebraically cancels mu3
% and makes m5 inert (a trap that cost a 45-minute sweep).
% n=701: amp and maxRe are converged there (deltas +0.72/+0.30/+1.37 dB, maxRe
% to 0.1). Map cleanliness is NOT, so the map columns are indicative only and
% any candidate is re-scored at n=1401.

M5 = [1 10 100 1000 1e4];
fprintf('\n  m=3 reference: amp d1 +81.15, d2 +84.17, maxRe +19.3\n');
fprintf('  m=4, chsz(4)=0.05 native, ohcgain=1, n=701\n\n');
fprintf('  sgn | m5 x  | amp d1  amp d2  amp d3 | maxRe     | range mono\n');
fprintf('%s\n', repmat('-',1,70));
R = struct('sgn',{},'m5',{},'S',{});
b5 = modpar26(4).m5o;
for sg = [1 -1]
    for k = 1:numel(M5)
        pa = modpar26(4);
        pa.ohcsgn = sg; pa.ohcgain = 1; pa.m5o = b5*M5(k);
        pa = setn(pa, 701);
        try
            S = score26(pa, 'fast', false);
        catch e
            fprintf('  %+3d | %5g | FAILED: %s\n', sg, M5(k), e.message); continue
        end
        fprintf('  %+3d | %5g | %+6.2f %+7.2f %+7.2f | %+9.1f | %5.2f %-4s\n', ...
                sg, M5(k), S.amp_gain, S.amp_d2, S.amp_d3, S.maxRe, ...
                S.bf_range, S.bf_mono);
        R(end+1).sgn=sg; R(end).m5=M5(k); R(end).S=S; %#ok<SAGROW>
        save('m4_bisect.mat','R');
    end
    fprintf('%s\n', repmat('-',1,70));
end
fprintf(['\n  READ: with sgn=-1, does amp rise toward m=3 (+81) as m5 freezes d3,\n' ...
         '  and does maxRe fall back into the healthy band (<= ~24)? If yes, the\n' ...
         '  d3 reaction path is the fault and it is localized. If amp stays low at\n' ...
         '  every m5, the inverted sign is not the whole story.\n']);
disp('M4_BISECT_DONE');
