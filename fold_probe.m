% WHERE IS THE 0.84 OCTAVE FOLD IN m=3b's MAP, AND IS IT PHYSICS OR DETECTION?
%
% CONTEXT. Fixing the apical edge artifact (find_bf returned the top bin at
% apical places, giving BF 24.98 kHz at levels of -39 dB) removed the
% fold==range signature but left a REAL fold: 0.84 oct in m=3b, 1.28 in the new
% m=4 default. That residual was masked all along and nothing in the project has
% examined it.
%
% HYPOTHESIS. It is the detector switching between two competing maxima in the
% d1 response, not a genuine map reversal. Resonance offsets measured earlier
% today: d1 is 28649 Hz at the base falling at -1.2683/cm; d2 is 14646 Hz
% falling at -1.7126/cm, so d2 sits 0.97 oct below d1 at the base and drifts to
% 3.21 oct below at the apex. If d1's response carries both the BM mode and the
% shear mode as separate peaks, the argmax will swap from one to the other where
% they cross, producing a step of roughly their separation. 0.84 oct is squarely
% in that range.
%
% THE TEST. Find the largest downward step, then at the places bracketing it
% report the TOP TWO local maxima of the same f-weighted curve local_tip uses.
%   CONFIRMED if the two peaks are present on BOTH sides and simply SWAP RANK,
%   with their separation matching the step size. The map is then fine and the
%   scalar BF is the wrong summary of a two-peaked response.
%   REFUTED if there is only one peak on either side, or the peak frequencies
%   themselves move discontinuously. That is a real reversal in the mechanics
%   and a genuine defect in the scaffolding model.
%
% METHOD copied from local_tip deliberately: same Am.*f weighting, same
% 0.15-18 kHz band, same interior XLO/XHI restriction, same edge rejection.
% Reinventing it would test a different detector than the one being diagnosed.

XLO = 0.05; XHI = 0.85;
ISVFRAC = (1391:-10:11)/1401;

pa = modpar26(3);
pa.isv = unique(max(1, min(pa.n, round(ISVFRAC*pa.n))), 'stable');
evalc('S = tdm26(0,pa,0,0);');

nf=numel(S.f); f=S.f(:); P=fft(S.ped); P=P(1:nf); np=size(S.d1,2);
xp=pa.isv(:)'/pa.n; inr=xp>=XLO & xp<=XHI;
bnd = f>0.15 & f<18; ib=find(bnd);
bf=nan(1,np); A=zeros(nf,np);
for i=find(inr)
    D=fft(S.d1(:,i)); a=abs(D(1:nf)./max(abs(P),eps)).*f; a(~bnd)=0;
    A(:,i)=a;
    [q,ip]=max(a); if (q<=0), continue; end
    if (ip<=ib(1) || ip>=ib(end)), continue; end     % edge rejection
    bf(i)=f(ip);
end

ok=find(isfinite(bf));
v=bf(ok); xv=xp(ok);
dd=diff(log2(v));
[worst,j]=min(dd);                                  % most negative = the fold
fprintf('\n  valid places %d of %d | range %.2f oct (%.2f - %.2f kHz)\n', ...
    numel(v), np, log2(max(v)/min(v)), min(v), max(v));
fprintf('  LARGEST DOWNWARD STEP %.3f oct, between x/L %.3f and %.3f\n', ...
    -worst, xv(j), xv(j+1));
fprintf('  BF %.3f kHz -> %.3f kHz (apex-to-base order)\n\n', v(j), v(j+1));

fprintf('  Top two local maxima of the SAME weighted curve at those places:\n');
fprintf('  x/L    | peak1 kHz  rel dB | peak2 kHz  rel dB | sep oct\n');
fprintf('%s\n', repmat('-',1,64));
for t = [j j+1]
    i = ok(t); a = A(:,i);
    lm = false(size(a)); lm(2:end-1) = a(2:end-1)>a(1:end-2) & a(2:end-1)>a(3:end);
    lm = lm & bnd;
    ii = find(lm);
    if (isempty(ii)), fprintf('  %.3f  | no interior local maximum\n', xp(i)); continue; end
    [~,srt] = sort(a(ii),'descend'); ii = ii(srt);
    f1 = f(ii(1)); a1 = a(ii(1));
    if (numel(ii)>=2)
        f2 = f(ii(2)); a2 = a(ii(2));
        fprintf('  %.3f  | %9.3f %7.2f | %9.3f %7.2f | %6.3f\n', ...
            xp(i), f1, 0.0, f2, 20*log10(max(a2,eps)/max(a1,eps)), ...
            abs(log2(f1/f2)));
    else
        fprintf('  %.3f  | %9.3f %7.2f | %9s %7s | %6s\n', ...
            xp(i), f1, 0.0, 'none', '-', '-');
    end
end
fprintf(['\n  CONFIRMED if both places show TWO peaks that swap rank and whose\n' ...
         '  separation matches the %.3f oct step: the map is sound and a scalar\n' ...
         '  BF is the wrong summary of a two-peaked response.\n' ...
         '  REFUTED if either place has a single peak, or the peak frequencies\n' ...
         '  themselves jump: that is a real reversal in the mechanics.\n'], -worst);
save('fold_probe.mat','bf','xp','A','f','ok','j');
disp('FOLD_PROBE_DONE');
