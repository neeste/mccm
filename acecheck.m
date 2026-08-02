% ACECHECK -- is the atlas's lone VIABLE verdict real, or metric quantization?
%
% THE CLAIM UNDER TEST. slopeatlas rated ace at RATE 4.12e-3 against a 3.4e-3
% threshold -- the only parameter of 16 to clear it, and sign-favorable (lowering
% the slope also lowers maxRe_osc). If true it is the lever for the long fit. Two
% reasons to distrust it before spending 12-20 h on it:
%
% 1. QUANTIZATION. aco and ace reported BIT-IDENTICAL dslope = -0.0685, i.e.
%    slope(+2%) - slope(-2%) = -0.00274 for both. Two unrelated parameters
%    landing on the same slope difference to four decimals is what a quantized
%    metric looks like, not matching physics. abr_metric's slope is a regression
%    on peak latencies picked off a dt=2e-6 grid, so a small perturbation can
%    leave every picked sample where it was. r3o, m1e and the exactly-0.0000
%    chsz pair are suspect for the same reason.
%
% 2. EXTRAPOLATION. ace is an EXPONENT: cp.ac = aco*exp(ace*x + acq*q), and it
%    sits at -0.399. The atlas rate says reaching the slope target needs a
%    normalized step of +2.82, which carries ace to +0.73 -- a SIGN FLIP of the
%    place-dependence of ac. A +-2% derivative extrapolated 282% is not evidence
%    of anything. The same arithmetic sends maperr to 90.97 - 331 < 0, which is
%    impossible, so the linearization is already known to break somewhere.
%
% So: sweep ace over the range the claim actually requires and look at the real
% numbers. Slope printed to 5 decimals so quantization is visible rather than
% inferred. Anchor row at the baseline value, as always.
t0=tic;
L=load('sweep_nch1b.mat'); pa0=L.R.pa;
vv=[pa0.ace -0.30 -0.20 -0.10 0 0.20 0.40 0.73];
fprintf('\n== ace sweep from sweep_nch1b (nch=1), cp.ac = aco*exp(ace*x) ==\n');
fprintf('   baseline ace %.5f: slope 0.60580  maperr 90.97  osc -56.6\n', pa0.ace);
fprintf('   atlas predicts ace=+0.73 reaches slope 0.413 AND improves osc.\n\n');
fprintf('     ace      slope     maperr   amp_gain   maxRe_osc  | verdict\n');
fprintf('   ----------------------------------------------------------------\n');
R=struct('ace',num2cell(vv));
for j=1:numel(vv)
    pa=pa0; pa.ace=vv(j);
    sl=NaN; mp=NaN; ag=NaN; mo=NaN;
    try, m=abr_metric(pa,false); if (m.ok), sl=m.slope; end, catch, end
    try
        S=score26(pa,'fast',false); mp=S.maperr; mo=S.maxRe_osc;
        if (isfield(S,'amp_gain')&&~isempty(S.amp_gain)), ag=double(S.amp_gain(1)); end
    catch, end
    v='';
    if (isfinite(mo) && mo>=0), v='UNSTABLE'; elseif (isfinite(mo) && mo>-40), v='inside margin'; end
    if (~isfinite(sl)), v=[v ' slope-NaN']; end
    fprintf('   %+7.4f  %8.5f  %8.2f   %+7.2f    %8.1f  | %s\n', vv(j), sl, mp, ag, mo, v);
    R(j).slope=sl; R(j).maperr=mp; R(j).amp=ag; R(j).osc=mo;
end
save('acecheck.mat','R','vv');
fprintf('\n   %.1f min.\n', toc(t0)/60);
fprintf('   REAL if slope falls toward 0.413 with osc staying below -40 and maperr sane.\n');
fprintf('   QUANTIZATION if slope barely moves across the whole range.\n');
disp('ACECHECK_DONE');
