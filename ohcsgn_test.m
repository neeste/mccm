% DOES ohcsgn=-1 MAKE THE 4-CHAMBER OHC PAIR ENERGY-INJECTING?
%
% The BM-work diagnostic showed, at the DEFAULT ohcsgn=+1:
%   ohcBM > 0 and scaling with gain  -> the force DOES reach and drive the BM
%   ohcP  < 0 at every stable gain   -> but the PAIR IS NET DISSIPATIVE
% Subtracting, the OC-height coordinate receives ohcP - ohcBM (negative): the
% force adds a little energy to the BM and removes more from the OC height, so
% it acts as a damper overall. That is why raising ohcgain buys a few dB and
% then instability instead of amplification. The fault is the SIGN/PHASE of the
% net energy flow, not the coupling path.
%
% HISTORY WORTH RESPECTING: an earlier session concluded from WNR MAGNITUDE that
% ohcsgn=-1 amplified, and the ohcP energy diagnostic contradicted it (showing
% dissipation / self-oscillation). WNR magnitude saturates and has misled this
% project three times, so it is NOT used as evidence here. The BM-work term did
% not exist for that earlier test, so this is a genuinely new measurement.
%
% PASS = ohcP > 0 (pair injects) AND ohcBM > 0 (BM specifically is driven) AND
% the model stays stable. ohcgain=0 is the validity control: both terms ~0.

b2 = modpar26(4).m2o; br1 = modpar26(4).r1o;
base = modpar26(4); base.m2o = b2*32; base.r1o = br1*0.5;

fprintf('\n  m=4, 2 kHz, 60 dB, m2o x32 + r1 x0.5\n');
fprintf('  sgn | ohcgain |      ohcP      |     ohcBM      |  lat ms | max|WNR|\n');
fprintf('%s\n', repmat('-',1,76));
for sg = [1 -1]
    for og = [0 1 2 4]
        p.fr = 2; p.lv = 60; p.pa = base;
        p.pa.ohcgain = og; p.pa.ohcsgn = sg; p.pa.hbmode = 'bm';
        try
            evalc('S = tdm26(''wnr1'', p, 0, 0);');
            d = S.dgn;
            oP = NaN; oB = NaN;
            if (isfield(d,'ohcP')),  oP = d.ohcP;  end
            if (isfield(d,'ohcBM')), oB = d.ohcBM; end
            fprintf('  %+3d | %7g | %+13.4e | %+13.4e | %6.2f | %9.3e\n', ...
                    sg, og, oP, oB, S.tpk, max(S.wnr));
        catch e
            fprintf('  %+3d | %7g | FAILED: %s\n', sg, og, e.message);
        end
    end
end

% ---- if sgn=-1 injects, score the amplifier and stability ------------------
fprintf('\n  Scoring both signs at ohcgain=1 (amp and maxRe, n=1401):\n');
fprintf('  sgn | amp dB  | maxRe     | range mono fold | contrast\n');
for sg = [1 -1]
    pa = base; pa.ohcgain = 1; pa.ohcsgn = sg;
    try
        S = score26(pa, 'fast', false);
        fprintf('  %+3d | %+6.2f | %+9.1f | %5.2f %-4s %.2f | %6.1f\n', ...
                sg, S.amp_gain, S.maxRe, S.bf_range, S.bf_mono, S.bf_fold, S.contrast);
    catch e
        fprintf('  %+3d | FAILED: %s\n', sg, e.message);
    end
end
fprintf(['\n  PASS = ohcP > 0 AND ohcBM > 0 AND stable. If sgn=-1 injects but the\n' ...
         '  model destabilizes, the pair is amplifying into self-oscillation and\n' ...
         '  the gain must be backed off rather than the sign reverted.\n']);
disp('OHCSGN_TEST_DONE');
