% CURRENT tip-tail status of the native 3-chamber, measured cleanly.
% The m=3 code path is unchanged this session (r4o=0 makes the r4 fix inert;
% m3form defaults to 0), so this is the model as it stands. m=2 for reference.
%
% tiptail_metric uses the d1-d2 differential on the ORIGINAL screened 7 places
% (0.6-16 kHz, CF-map screened, R2-gated) -- the same metric used for every
% prior d number, so these are directly comparable to the 0.784 / 0.445 history.
% Also reports the native gam and whether the operating point is the degenerate
% one flagged in the ohc_tip sweep (BM d1 apical blow-up at gam=1).

for nch=[2 3]
    pa=modpar26(nch);
    fprintf('\n===== m=%d  (native gam=%.3f) =====\n', nch, pa.gam);
    m=tiptail_metric(pa,false);
    if (m.ok || m.nvalid>=4)
        fprintf('  d = %.3f   R2 = %.2f   nvalid = %d\n', m.d, m.r2, m.nvalid);
        fprintf('  BF (kHz)   : %s\n', num2str(m.BF,'%7.2f'));
        fprintf('  contrast dB: %s   (tip-tail per place)\n', num2str(m.contrast,'%7.1f'));
        fprintf('  chi (hi-CF tip-tail) = %.1f dB\n', m.chi);
    else
        fprintf('  tiptail_metric NOT OK: %s  (nvalid=%d)\n', m.msg, m.nvalid);
    end
    if (~isempty(m.msg)), fprintf('  note: %s\n', m.msg); end
end

% Cross-check the BM-observable operating point at native gam (the ohc_tip view):
% is the native 3-chamber sitting on the degenerate apical resonance?
fprintf('\n----- BM-observable check, m=3 native gam -----\n');
pa=modpar26(3); ISV=1391:-10:11; XLO=0.05; XHI=0.85;
evalc('S=tdm26(0,pa,0,0);');
nf=numel(S.f); f=S.f(:); P=fft(S.ped); P=P(1:nf); np=size(S.d1,2);
xp=pa.isv(:)'/pa.n; inr=xp>=XLO & xp<=XHI;
[~,ifq]=min(abs(f-2)); A2=zeros(1,np);
for i=1:np, D=fft(S.d1(:,i)); A2(i)=abs(D(ifq))/max(abs(P(ifq)),eps); end
[pa_all,ia]=max(A2); A2i=A2; A2i(~inr)=0; [pa_in,ii]=max(A2i);
locw={'interior','BOUNDARY/apical'};
fprintf('  |H_BM|@2k: peak_all at x/L=%.3f (%s)   peak_interior at x/L=%.3f\n', ...
        xp(ia), locw{(~inr(ia))+1}, xp(ii));
fprintf('  (native gam=%.2f; if the peak is apical/boundary, the BM view is the\n', pa.gam);
fprintf('   degenerate resonance even though the d1-d2 metric reads a clean d.)\n');
disp('TT_STATUS_DONE');
