% DOES THE OHC FORCE REACH THE BM IN THE 4-CHAMBER?
%
% ohcP  = -sum(act.*(vc-v1))  work on the RELATIVE (BM-to-RL) coordinate. This
%         is the correct measure for an internal force pair, and it is what the
%         existing energy diagnostic reports. It does NOT say whether the BM is
%         amplified.
% ohcBM = sum(act.*v1)        work against BM VELOCITY. The BM receives +act as
%         the reaction in s1, so this is the power delivered to BM motion. This
%         is the quantity that decides amplification, and it was never measured.
%
% READING:
%   ohcBM > 0 and scaling with ohcgain  -> the force does drive the BM; the low
%       gain has some other cause (mass ratio, phase, saturation).
%   ohcBM ~ 0 while ohcP > 0            -> the force injects energy into the
%       OC-height coordinate only and never reaches the BM. That explains an
%       amplifier stuck near 2.5 dB regardless of gain or damping, and points at
%       the DRIVING COORDINATE (act is driven by d2, the SV-SS shear, whereas
%       m=1/m=2 drive the active force from d1-d2, which contains BM motion).
%   ohcBM < 0                           -> the force DAMPS the BM.
%
% The wnr1 path is used because dgn (which carries ohcP/ohcBM) is built there.

b2 = modpar26(4).m2o; br1 = modpar26(4).r1o;
base = modpar26(4); base.m2o = b2*32; base.r1o = br1*0.5;

fprintf('\n  m=4, 2 kHz, 60 dB, m2o x32 + r1 x0.5\n');
fprintf('  ohcgain |      ohcP      |     ohcBM     |  ratio BM/rel | lat ms | max|WNR|\n');
fprintf('%s\n', repmat('-',1,80));
for og = [0 1 2 4 8]
    p.fr = 2; p.lv = 60; p.pa = base; p.pa.ohcgain = og; p.pa.hbmode = 'bm';
    try
        evalc('S = tdm26(''wnr1'', p, 0, 0);');
        d = S.dgn;
        oP  = NaN; oB = NaN;
        if (isfield(d,'ohcP')),  oP = d.ohcP;  end
        if (isfield(d,'ohcBM')), oB = d.ohcBM; end
        rat = NaN; if (isfinite(oP) && abs(oP) > 0), rat = oB/oP; end
        fprintf('  %7g | %+13.4e | %+13.4e | %13.4g | %6.2f | %9.3e\n', ...
                og, oP, oB, rat, S.tpk, max(S.wnr));
    catch e
        fprintf('  %7g | FAILED: %s\n', og, e.message);
    end
end
fprintf(['\n  ohcgain=0 is the control: both terms must be ~0 there (act is scaled\n' ...
         '  to zero), which validates the measurement before the other rows are read.\n']);
disp('BM_WORK_DIAG_DONE');
