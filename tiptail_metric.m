function m = tiptail_metric(pa, verbose)
%TIPTAIL_METRIC  Forward-latency exponent d and tip-tail contrast from a
%   LINEARIZED click (short-chirp) response.  ~11 s for m=2, versus ~100 s for
%   the tone-burst metric, so it supports a much larger fit budget.
%
%   Delay accumulates in the TIP, and d<1 requires delay-in-cycles N=tau*BF to
%   RISE with CF.  m.d is fitted from log N vs log BF.
%
%   CF-MAP SCREENING is mandatory: a place whose BF breaks monotonicity is a
%   map failure, not tuning, and silently wrecks the regression (it produced a
%   spurious d=0.183 for the native 3-chamber).  isv runs apical->basal, so BF
%   must increase along it.
if (nargin<2), verbose=false; end
m.ok=false; m.msg=''; m.d=NaN; m.r2=NaN; m.BF=[]; m.N=[]; m.contrast=[]; m.Q3=[]; m.nvalid=0; m.chi=NaN;
try
    if (verbose), S=tdm26(0,pa,0,0); else, evalc('S=tdm26(0,pa,0,0);'); end
catch e
    m.msg=sprintf('tdm26 threw: %s',e.message); return
end
nf=numel(S.f); f=S.f(:); P=fft(S.ped); P=P(1:nf);
np=size(S.d1,2); BF=nan(1,np); CO=nan(1,np); Q3=nan(1,np); NN=nan(1,np);
for i=1:np
    D=fft(S.d1(:,i)-S.d2(:,i)); D=D(1:nf);
    H=D./max(abs(P),eps); A=abs(H).*f(:);
    A(~(f>0.15 & f<18))=0;
    [pk,ip]=max(A); if (pk<=0), continue; end
    bf=f(ip); if (bf<0.6 || bf>16), continue; end   % low-CF phase unwrap is
                                                     % unreliable in a click run
    [~,it]=min(abs(f-bf/4));
    CO(i)=20*log10(pk/max(A(it),eps));
    Q3(i)=bf/max(deal_bw(A,f,ip,pk/10^(3/20)),eps);
    ph=unwrap(angle(H)); w=2*pi*f*1000; gd=-gradient(ph,w);
    NN(i)=median(gd(max(1,ip-3):min(nf,ip+3)))*1000*bf;
    BF(i)=bf;
end
% --- CF-map screening: keep the longest run of strictly increasing BF ---
ok=isfinite(BF)&isfinite(NN)&NN>0;
idx=find(ok); best=[]; run=[];
for k=1:numel(idx)
    if (isempty(run) || BF(idx(k))>BF(run(end))*1.05), run(end+1)=idx(k); %#ok<AGROW>
    else, if (numel(run)>numel(best)), best=run; end; run=idx(k); end
end
if (numel(run)>numel(best)), best=run; end
m.nvalid=numel(best);
if (m.nvalid<4), m.msg=sprintf('only %d valid places after CF-map screening',m.nvalid); return; end
m.BF=BF(best); m.N=NN(best); m.contrast=CO(best); m.Q3=Q3(best);
b=polyfit(log(m.BF),log(m.N),1);
m.d=1-b(1);
r=log(m.N)-polyval(b,log(m.BF));                 % fit quality: a bad place shows
ss=sum((log(m.N)-mean(log(m.N))).^2);            % up as low R2, so it cannot
m.r2=1-sum(r.^2)/max(ss,eps);                    % corrupt d unnoticed
if (m.r2<0.80), m.msg=sprintf('poor N-vs-BF fit (R2=%.2f)',m.r2); end
m.chi=m.contrast(end);          % tip-tail contrast at the HIGHEST valid CF
m.ok=true;
end
