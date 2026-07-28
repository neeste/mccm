% BM-RIGIDITY, done properly: KINEMATIC clamp (pa.bmrigid), not stiffening.
% Stiffening k1o broke the explicit time step and returned NaN for every chamber
% count, measuring the integrator rather than the physics.
% CLAIM UNDER TEST: m=1/m=2 micromechanics vanish with a rigid BM (the BM is the
% only partition the fluid acts on); m>=4 should survive because the drive
% reaches d2 through the SV-SS pressure difference without needing d1 to move.
ISV=1391:-10:11; XLO=0.05; XHI=0.85;
fprintf('\n  m | BM     | pk|d1|    pk|d2|    pk|d3|   | amp d2  amp d3\n');
fprintf('%s\n',repmat('-',1,66));
for m=[1 2 3 4]
  for rg=[0 1]
    pa=modpar26(m); pa.bmrigid=rg; pa.isv=ISV;
    [pk,g]=probe2(pa,ISV,XLO,XHI);
    lbl={'free  ','RIGID '};
    fprintf('  %d | %s | %9.2e %9.2e %9.2e | %+6.2f %+7.2f\n', ...
            m, lbl{rg+1}, pk(1), pk(2), pk(3), g(2), g(3));
  end
  fprintf('%s\n',repmat('-',1,66));
end
fprintf('\n  RIGID rows: does d2/d3 still MOVE and still AMPLIFY?\n');
disp('BM_RIGID2_DONE');

function [pk,gain]=probe2(pa,ISV,XLO,XHI)
pk=nan(1,3); gain=nan(1,3);
p1=pa; p0=pa;
if (pa.m>=4), og=1; if(isfield(pa,'ohcgain')),og=pa.ohcgain;end
    p1.ohcgain=og; p0.ohcgain=0;
else, p1.gam=pa.gam; p0.gam=0; end
A=amps2(p1,ISV,XLO,XHI); B=amps2(p0,ISV,XLO,XHI);
pk=A.pkabs;
for q=1:3
    if(isfinite(A.pk2k(q))&&isfinite(B.pk2k(q))&&B.pk2k(q)>0)
        gain(q)=20*log10(A.pk2k(q)/B.pk2k(q));
    end
end
end
function R=amps2(pa,ISV,XLO,XHI)
R.pkabs=nan(1,3); R.pk2k=nan(1,3); pa.isv=ISV;
try, evalc('S=tdm26(0,pa,0,0);'); catch, return; end
if(any(~isfinite(S.d1(:)))), return; end
nf=numel(S.f); f=S.f(:); P=fft(S.ped); P=P(1:nf); np=size(S.d1,2);
xp=pa.isv(:)'/pa.n; inr=xp>=XLO&xp<=XHI;
D={S.d1,S.d2,[]};
if(isfield(S,'d3')&&~isempty(S.d3)&&any(S.d3(:)~=0)), D{3}=S.d3; end
[~,jf]=min(abs(f-2));
for q=1:3
    if(isempty(D{q})), continue; end
    R.pkabs(q)=max(abs(D{q}(:)));
    Am=zeros(1,np);
    for i=1:np, Z=fft(D{q}(:,i)); Am(i)=abs(Z(jf)/max(abs(P(jf)),eps)); end
    Am(~isfinite(Am))=0; Am(~inr)=0;
    if(max(Am)>0), R.pk2k(q)=max(Am); end
end
end
