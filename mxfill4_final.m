% VALIDATE THE m=4 mxfill EXTENSION -- using the new in-file cfmap entry point.
%
% Four external harnesses failed because fdm26's internals (fdmod23, tuning_err,
% imped, mxfill) are LOCAL subfunctions. fdm26(struct('cfmap',1,...)) now returns
% the peak PLACE per frequency, computed the way tuning_err does it. Comparison:
%   fdm26 peak place (Db, BM)  vs  tdm26 peak place (d1 click, BM)
% over the maperr frequency grid. m=2 and m=3 CALIBRATE the comparison; accept
% m=4 only if comparable. Two sanity gates so a fifth harness bug announces
% itself instead of hiding in an rms:
%   (a) fdm26 peak place must move BASALLY (x/L decreasing) with frequency;
%   (b) tdm26 must do the same, scored on the INTERIOR only (x/L 0.05..0.85) --
%       the dense grid otherwise pins to the helicotrema/stapes boundaries.

flmap = 500*2.^(-1:5);              % 250..16000 Hz (the maperr grid)
ISV   = 1391:-10:11; XLO=0.05; XHI=0.85;

fprintf('\n  nch T4sgn |  peak place x/L : fdm26 (Db) over tdm26 (d1)        | rms  | gates\n');
fprintf('%s\n', repmat('-',1,100));
for nch=[2 3 4]
  sgl=1; if (nch==4), sgl=[1 -1]; end
  for sg=sgl
    pa=modpar26(nch); if (nch==4), pa.T4sgn=sg; end
    xf=nan(1,numel(flmap)); xt=xf;
    try
        Rc=fdm26(struct('cfmap',1,'pa',pa,'flst',flmap));
        xf=Rc.xpk_bm;
    catch e
        fprintf('  %d  %+4d | cfmap FAILED: %s\n', nch, sg, e.message); continue
    end
    try
        pt=pa; pt.isv=ISV;
        evalc('S=tdm26(0,pt,0,0);');
        nf=numel(S.f); f=S.f(:); P=fft(S.ped); P=P(1:nf); np=size(S.d1,2);
        xp=pt.isv(:)'/pt.n; inr = xp>=XLO & xp<=XHI;
        H=zeros(nf,np);
        for i=1:np, D=fft(S.d1(:,i)); H(:,i)=D(1:nf)./max(abs(P),eps); end
        Am=abs(H);
        for k=1:numel(flmap)
            [~,jf]=min(abs(f-flmap(k)/1000));
            e2=Am(jf,:); e2(~isfinite(e2))=0; e2(~inr)=0;
            [pk,jb]=max(e2); if (pk>0), xt(k)=xp(jb); end
        end
    catch e
        fprintf('  %d  %+4d | tdm26 FAILED: %s\n', nch, sg, e.message); continue
    end
    d=xf-xt; rms=sqrt(mean(d(isfinite(d)).^2));
    g=''; vf=xf(isfinite(xf)); vt=xt(isfinite(xt));
    if (numel(vf)<3 || ~all(diff(vf)<=0)), g=[g ' FDM-NONMONO']; end
    if (numel(vt)<3 || ~all(diff(vt)<=0)), g=[g ' TDM-NONMONO']; end
    if (isempty(g)), g=' ok'; end
    fprintf('  %d  %+4d | f: %s\n', nch, sg, num2str(xf,'%7.3f'));
    fprintf('            | t: %s | %.4f |%s\n', num2str(xt,'%7.3f'), rms, g);
  end
end
fprintf('\n  freqs (Hz): %s\n', num2str(flmap,'%7.0f'));
fprintf(['\nACCEPT m=4 only if (i) both gates read ok for m=2 and m=3, establishing\n' ...
         'the comparison is sound, and (ii) the m=4 rms is comparable to theirs.\n' ...
         'The smaller-rms T4sgn is the correct sign.\n']);
disp('MXFILL4_FINAL_DONE');
