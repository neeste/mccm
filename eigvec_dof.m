% WHICH DOF DOES THE m=4 STATIC DIVERGENCE OCCUPY?
% coupeig (with eigvec=1) now returns the max-real eigenvector, decomposed into
% BM / shear / OC-height participation and each DOF's peak place. The DOF with
% the dominant fraction is the divergent coordinate -- i.e. the one that ACTUALLY
% needs the restoring force, measured rather than guessed. (The k5 sweep showed
% stiffening OC-height made it worse, so the prediction is the mode does NOT live
% mainly in OC-height.)
%
% Configs: the strongly-unstable gam=0.5 (maxRe +16862) and gam=0.3, plus the
% stable native gam=1.0 as a control -- its max-real mode should be the spurious
% +0.0 boundary mode, not a physical instability, and should look different.

cfgs = { 0.50, 'gam=0.50  (collapse, maxRe~+16862)'
         0.30, 'gam=0.30  (diverges)'
         1.00, 'gam=1.00  (stable native -- control)' };

fprintf('\n  config                              maxRe(all)  |  BM     shear   OChgt  |  ME   | peak x/L (BM/sh/OC)\n');
fprintf('%s\n', repmat('-',1,108));
for c=1:size(cfgs,1)
    gam=cfgs{c,1}; lbl=cfgs{c,2};
    pa=modpar26(4); pa.gam=gam;
    try
        evalc('E=tdm26(''coupeig'',struct(''pa'',pa,''eigvec'',1));');
        f=E.uDOFfrac; x=E.uDOFpkx;
        fprintf('  %-34s  %+10.1f  | %5.2f  %5.2f  %5.2f  | %5.2f | %.2f / %.2f / %.2f\n', ...
                lbl, E.maxRe, f(1), f(2), f(3), E.uMEfrac, x(1), x(2), x(3));
        fprintf('       top eigenvalue: Re=%+.1f  f=%.3f kHz\n', real(E.ulam), abs(imag(E.ulam))/2/pi/1000);
    catch e
        fprintf('  %-34s  FAILED: %s\n', lbl, e.message);
    end
end
fprintf('\nDominant fraction = the divergent DOF. peak x/L: higher = more apical.\n');
fprintf('If BM or shear dominates (not OChgt), that explains why stiffening k5\n');
fprintf('(OC-height) did not help -- the restoring force belongs on another coord.\n');
disp('EIGVEC_DOF_DONE');
