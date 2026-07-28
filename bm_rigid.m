% DOES AMPLIFICATION SURVIVE A RIGID BASILAR MEMBRANE?  (SN, 2026-07-26)
%
% Lizards, birds and amphibians have no basilar membrane in the mammalian sense
% yet show sharp tuning and active emissions, so amplification is a property of
% the hair cell and tectorial micromechanics rather than of the BM. A model that
% captures that should still amplify when the BM is made rigid.
%
% THE STRUCTURAL CLAIM UNDER TEST. In m=1/m=2 the BM is the ONLY partition the
% fluid acts on, so clamping it removes the drive to everything downstream and
% the micromechanics go quiet. In m>=4 the drive reaches d2 through the SV-SS
% pressure difference, which does not require d1 to move. If d2/d3 survive
% d1 -> 0 in m=4 but every coordinate dies in m=1/m=2, the 4-chamber reproduces
% the BM-free amplification mode and the 1-chamber structurally cannot. That
% would be a clear qualitative advantage, and unlike "better fit" it is a claim
% the 4-chamber can actually support.
%
% METHOD. Stiffen the BM by scaling k1o. This pushes the CF map far up, which is
% expected and irrelevant here: the question is whether the OTHER coordinates
% still respond and still amplify, not whether the tonotopy is right.
% Reported per coordinate: peak |d| (does it move at all) and active-vs-passive
% gain (does the amplifier still work). d3 exists only for m>=4.
%
% NOTE d2 and d3 are ABSOLUTE displacements: a(i2)=a(i1)+s2/m2 and
% a(i3)=a(i1)+s3/m5 both add the BM acceleration. So d3 is the absolute
% RETICULAR LAMINA position, which is the coordinate the mouse OCT work says may
% matter more than the BM.

KS = [1 1e2 1e4 1e6];
MS = [1 2 3 4];
ISV = 1391:-10:11; XLO=0.05; XHI=0.85;

fprintf('\n  BM stiffened by k1o x KS. peak|d| and active-vs-passive gain per DOF.\n');
fprintf('  m | k1o x | pk|d1|    pk|d2|    pk|d3|   | amp d1  amp d2  amp d3\n');
fprintf('%s\n', repmat('-',1,78));
for m = MS
    base = modpar26(m);
    for ks = KS
        pa = base; pa.k1o = base.k1o*ks; pa.isv = ISV;
        [p1,g1] = probe(pa, ISV, XLO, XHI);
        fprintf('  %d | %5.0e | %9.2e %9.2e %9.2e | %+6.2f %+7.2f %+7.2f\n', ...
                m, ks, p1(1), p1(2), p1(3), g1(1), g1(2), g1(3));
    end
    fprintf('%s\n', repmat('-',1,78));
end
fprintf(['\n  READ: if m=1/m=2 lose ALL coordinates as k1o rises while m=4 keeps\n' ...
         '  d2/d3 motion AND d2/d3 gain, the 4-chamber supports BM-free\n' ...
         '  amplification. If m=4 also goes quiet, its micromechanics are not\n' ...
         '  independent of the BM and that advantage claim is not available.\n']);
disp('BM_RIGID_DONE');

function [pk, gain] = probe(pa, ISV, XLO, XHI)
pk = nan(1,3); gain = nan(1,3);
pa1 = pa; pa0 = pa;
if (pa.m >= 4)
    og = 1; if (isfield(pa,'ohcgain')), og = pa.ohcgain; end
    pa1.ohcgain = og; pa0.ohcgain = 0;
else
    pa1.gam = pa.gam; pa0.gam = 0;
end
A = amps(pa1, ISV, XLO, XHI);
B = amps(pa0, ISV, XLO, XHI);
pk = A.pkabs;
for q = 1:3
    if (isfinite(A.pk2k(q)) && isfinite(B.pk2k(q)) && B.pk2k(q) > 0)
        gain(q) = 20*log10(A.pk2k(q)/B.pk2k(q));
    end
end
end

function R = amps(pa, ISV, XLO, XHI)
R.pkabs = nan(1,3); R.pk2k = nan(1,3);
pa.isv = ISV;
try, evalc('S=tdm26(0,pa,0,0);'); catch, return; end
if (any(~isfinite(S.d1(:)))), return; end
nf=numel(S.f); f=S.f(:); P=fft(S.ped); P=P(1:nf); np=size(S.d1,2);
xp=pa.isv(:)'/pa.n; inr = xp>=XLO & xp<=XHI;
D = {S.d1, S.d2, []};
if (isfield(S,'d3') && ~isempty(S.d3) && any(S.d3(:)~=0)), D{3} = S.d3; end
[~,jf] = min(abs(f-2));
for q = 1:3
    if (isempty(D{q})), continue; end
    R.pkabs(q) = max(abs(D{q}(:)));
    Am = zeros(1,np);
    for i = 1:np
        Z = fft(D{q}(:,i)); Am(i) = abs(Z(jf)/max(abs(P(jf)),eps));
    end
    Am(~isfinite(Am)) = 0; Am(~inr) = 0;
    if (max(Am) > 0), R.pk2k(q) = max(Am); end
end
end
