% WHY DOES m3form=1 RETURN NaN AT LOW gam (gam=0, 0.3)?
% Two very different causes, separated here:
%   (I)  DIVERGENCE  -- the time march blew up; d1 contains Inf/NaN.
%   (II) NO-PEAK     -- the run is finite but the excitation pattern has no
%                       interior CF peak (peak pinned to a boundary, or no
%                       traveling-wave structure at all).
% Also runs coupeig for the linear-stability verdict, and compares against the
% d2-only form at the SAME gam so the diagnosis is differential.
%
% gam=0 is the crucial row: gam multiplies s4tmp, so gam=0 is the PURELY PASSIVE
% model. If the passive baseline is already broken, the active form is moot.

ISV = 1391:-10:11; XLO=0.05; XHI=0.85;
for gam = [0 0.30 0.70]
for mf = [0 1]
    pa = modpar26(3); pa.gam = gam; pa.isv = ISV; pa.m3form = mf;
    fprintf('\n===== m=3  gam=%.2f  m3form=%d =====\n', gam, mf);
    % ---- run and check finiteness -------------------------------------------
    err='';
    try, evalc('S = tdm26(0,pa,0,0);'); catch e, err=e.message; end
    if (~isempty(err)), fprintf('  tdm26 THREW: %s\n', err); continue; end
    d1=S.d1; nfin=sum(~isfinite(d1(:)));
    fprintf('  d1: %d/%d non-finite   max|d1|=%.3e   ped max=%.3e\n', ...
            nfin, numel(d1), max(abs(d1(isfinite(d1)))), max(abs(S.ped)));
    if (nfin>0)
        [rr,~]=find(~isfinite(d1),1); fprintf('  --> DIVERGENCE (first non-finite at time-row %d)\n', rr);
    end
    % ---- excitation pattern at 2 kHz: where is the peak? --------------------
    nf=numel(S.f); f=S.f(:); P=fft(S.ped); P=P(1:nf);
    np=size(d1,2); xp=pa.isv(:)'/pa.n; inr = xp>=XLO & xp<=XHI;
    [~,ifq]=min(abs(f-2)); H2=zeros(1,np);
    for i=1:np, D=fft(d1(:,i)); H2(i)=D(ifq)/max(abs(P(ifq)),eps); end
    A2=abs(H2);
    [pkall,iall]=max(A2);                 % peak over ALL places
    A2i=A2; A2i(~inr)=0; [pkin,iin]=max(A2i);   % peak in interior
    locw={'interior','BOUNDARY'};
    fprintf('  |H| @2k: peak_all=%.3e at x/L=%.3f (%s)   peak_interior=%.3e at x/L=%.3f\n', ...
            pkall, xp(iall), locw{(~inr(iall))+1}, pkin, xp(iin));
    % how many interior places have a real BF peak (the CF-map count)
    nbf=0;
    for i=find(inr)
        a=zeros(nf,1);
        for kk=1:1, end
        D=fft(d1(:,i)); Hc=D(1:nf)./max(abs(P),eps); aa=abs(Hc).*f; aa(~(f>0.15&f<18))=0;
        [pk,ip]=max(aa); if (pk>0 && f(ip)>=0.5 && f(ip)<=16.5), nbf=nbf+1; end
    end
    fprintf('  interior places with a valid BF peak: %d / %d\n', nbf, sum(inr));
    % ---- linear stability ---------------------------------------------------
    try
        evalc('E=tdm26(''coupeig'',struct(''pa'',pa));');
        v='sub-critical (stable)'; if (E.maxRe_osc>=0), v='UNSTABLE'; end
        fprintf('  coupeig maxRe_osc = %+.2f  (%s)\n', E.maxRe_osc, v);
    catch e
        fprintf('  coupeig failed: %s\n', e.message);
    end
end
end
disp('M3FORM_NAN_DONE');
