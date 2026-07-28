% IS THE TRACKER FOLLOWING THE SUB-DOMINANT MODE?
%
% SYMPTOM. With the continuity tracker, contrast went NEGATIVE for both m=4
% configs (-13.6, -14.4) while m=3b held at 9.2, and bf_hi doubled from 7.40 to
% 14.77 kHz. Negative contrast means the response at bf/4 EXCEEDS the response
% at bf, i.e. the reported peak is not the largest one.
%
% THE LIKELY CAUSE, and it is conceptual rather than a coding slip. Continuity
% and CF are different quantities. At a mode crossing the mode trajectory
% continues smoothly while the DOMINANT peak switches branches. A continuity
% tracker therefore follows the weaker mode past the crossing by design. The CF
% map is defined by the largest response, so the argmax jump was never a
% detector failure -- it was the model reporting two modes crossing with
% near-equal amplitude (fold_probe measured the margin at 0.01 dB).
%
% THE TEST. At every valid place, compare the TRACKED peak's amplitude against
% the LARGEST peak's amplitude in the same weighted curve. If the tracked peak
% is materially smaller over a contiguous basal stretch, the tracker has left
% the dominant mode and the continuity approach is answering the wrong question.
%
% This distinguishes two repairs:
%   tracker is simply buggy      -> a few scattered places disagree
%   tracker is conceptually wrong-> a contiguous run disagrees, starting at the
%                                   crossing and continuing to the array end

XLO=0.05; XHI=0.85; ISVFRAC=(1391:-10:11)/1401;
pa = modpar26(4);
pa.isv = unique(max(1,min(pa.n,round(ISVFRAC*pa.n))),'stable');
evalc('S = tdm26(0,pa,0,0);');

nf=numel(S.f); f=S.f(:); P=fft(S.ped); P=P(1:nf); np=size(S.d1,2);
xp=pa.isv(:)'/pa.n; inr=xp>=XLO & xp<=XHI;
bnd=f>0.15&f<18; ib=find(bnd);
pkf=cell(1,np); pka=cell(1,np);
for i=find(inr)
    a=abs(fft(S.d1(:,i))); a=a(1:nf)./max(abs(P),eps).*f; a(~bnd)=0;
    lm=false(nf,1); lm(2:nf-1)=a(2:nf-1)>a(1:nf-2)&a(2:nf-1)>a(3:nf);
    lm=lm&bnd; lm(ib(1))=false; lm(ib(end))=false;
    jj=find(lm); if(isempty(jj)), continue; end
    [~,sr]=sort(a(jj),'descend'); jj=jj(sr);
    pkf{i}=f(jj); pka{i}=a(jj);
end
have=find(~cellfun(@isempty,pkf));
mrg=-inf(1,np);
for i=have
    p=pka{i}; if(numel(p)>=2), mrg(i)=20*log10(p(1)/max(p(2),eps)); else, mrg(i)=60; end
end
[~,i0]=max(mrg);
bf=nan(1,np); rank1=nan(1,np); deficit=nan(1,np);
bf(i0)=pkf{i0}(1); rank1(i0)=1; deficit(i0)=0;
for dir=[1 -1]
    ref=bf(i0); i=i0+dir;
    while (i>=1 && i<=np)
        if (~isempty(pkf{i}))
            [dl,kk]=min(abs(log2(pkf{i}/ref)));
            if (dl<=0.5)
                bf(i)=pkf{i}(kk); ref=bf(i); rank1(i)=kk;
                deficit(i)=20*log10(pka{i}(kk)/pka{i}(1));  % 0 if tracking the largest
            end
        end
        i=i+dir;
    end
end
ok=find(isfinite(bf));
fprintf('\n  seed at x/L %.3f (margin %.2f dB over runner-up)\n', xp(i0), mrg(i0));
fprintf('  valid places %d | tracking the LARGEST peak at %d of them (%.0f%%)\n', ...
    numel(ok), sum(rank1(ok)==1), 100*mean(rank1(ok)==1));
bad = ok(rank1(ok)~=1);
if (isempty(bad))
    fprintf('  tracker never leaves the dominant mode -- symptom is elsewhere\n');
else
    fprintf('  leaves the dominant mode at %d places, x/L %.3f to %.3f\n', ...
        numel(bad), min(xp(bad)), max(xp(bad)));
    fprintf('  worst amplitude deficit %.2f dB\n', min(deficit(bad)));
    d=diff(sort(bad)); runs=1+sum(d>1);
    fprintf('  those places form %d contiguous run(s)\n', runs);
end
fprintf('\n  x/L    | tracked kHz | rank | deficit dB\n');
fprintf('%s\n', repmat('-',1,48));
sel = ok(max(1,numel(ok)-11):end);      % the most BASAL dozen
for i = sel
    fprintf('  %.3f  | %11.3f | %4d | %10.2f\n', xp(i), bf(i), rank1(i), deficit(i));
end
fprintf(['\n  A CONTIGUOUS basal run with rank>1 means continuity is answering the\n' ...
         '  wrong question and the argmax (plus the edge fix) is the correct CF\n' ...
         '  definition, with the fold reporting genuine mode degeneracy.\n' ...
         '  SCATTERED disagreements would instead mean the tracker is merely buggy.\n']);
disp('TRACKER_DIAG_DONE');
