% VALIDATE THE m=4 mxfill EXTENSION -- third attempt, using fdm26's OWN path.
%
% The two previous harnesses guessed at fdm26's output conventions and both
% failed identically for m=2/m=3, proving the harness wrong rather than the
% physics. This one uses fdmod23 (fdm26.m:522), the function fdm26's own
% tuning_err calls: it returns Db/Dh as (PLACE x FREQUENCY), and tuning_err
% takes max over PLACE to locate the CF place. So the comparison is
%   PEAK PLACE vs FREQUENCY -- fdm26 (Db) against tdm26 (d1 click) --
% both on BM displacement, which is apples-to-apples.
%
% m=2 and m=3 are the CALIBRATION: fdm26 is long-established there, so their
% agreement sets the bar. Accept m=4 only if comparable. If m=2/m=3 disagree
% badly, the harness is STILL wrong and the m=4 number means nothing.

flmap = 500*2.^(-1:5);          % 250..16000 Hz, the same grid maperr uses
ISV   = 1391:-10:11;            % dense places for tdm26

fprintf('\n  nch T4sgn |  peak place x/L per freq (fdm26 Db / tdm26 d1)      | rms dx/L\n');
fprintf('%s\n', repmat('-',1,94));
for nch=[2 3 4]
  sgl=1; if (nch==4), sgl=[1 -1]; end
  for sg=sgl
    pa=modpar26(nch); if (nch==4), pa.T4sgn=sg; end
    xf=nan(1,numel(flmap)); xt=xf;
    % ---- fdm26 peak place via the fdsolve dispatch ----
    % fdmod23/tuning_err/imped/mxfill are LOCAL subfunctions of fdm26.m and are
    % NOT callable from here (that was the previous failure). fd_solve is
    % reachable through the 'fdsolve' dispatch and returns R.Db, the same BM
    % displacement fdmod23 provides. Reduce by max over PLACE (as tuning_err
    % does), NOT over frequency -- reducing over frequency was failure #1.
    try
        for k=1:numel(flmap)
            Rj=fdm26(struct('fdsolve',1,'pa',pa,'f',flmap(k)));
            if (~isstruct(Rj) || ~isfield(Rj,'Db')), error('no Db from fdsolve'); end
            a=abs(Rj.Db(:,1)); a(~isfinite(a))=0;
            [pk,ip]=max(a); if (pk>0), xf(k)=ip/pa.n; end
        end
    catch e
        fprintf('  %d  %+4d | fdsolve FAILED: %s\n', nch, sg, e.message); continue
    end
    % ---- tdm26 peak place from the linear click ----
    try
        pt=pa; pt.isv=ISV;
        evalc('S=tdm26(0,pt,0,0);');
        nf=numel(S.f); f=S.f(:); P=fft(S.ped); P=P(1:nf); np=size(S.d1,2);
        xp=pt.isv(:)'/pt.n;
        H=zeros(nf,np);
        for i=1:np, D=fft(S.d1(:,i)); H(:,i)=D(1:nf)./max(abs(P),eps); end
        Am=abs(H);
        for k=1:numel(flmap)
            [~,jf]=min(abs(f-flmap(k)/1000));
            e2=Am(jf,:); e2(~isfinite(e2))=0;
            [pk,jb]=max(e2); if (pk>0), xt(k)=xp(jb); end
        end
    catch e
        fprintf('  %d  %+4d | tdm26 FAILED: %s\n', nch, sg, e.message); continue
    end
    d=xf-xt; rms=sqrt(mean(d(isfinite(d)).^2));
    % SANITY: peak place must move BASALLY (x/L smaller) as frequency rises.
    % If it does not, the harness is wrong again and the rms is meaningless.
    mono=''; v=xf(isfinite(xf));
    if (numel(v)>=3 && ~all(diff(v)<0)), mono='  <-- NOT MONOTONIC (harness suspect)'; end
    fprintf('  %d  %+4d | f: %s%s\n', nch, sg, num2str(xf,'%7.3f'), mono);
    fprintf('            | t: %s |  %.4f\n', num2str(xt,'%7.3f'), rms);
  end
end
fprintf(['\n  freqs (Hz): %s\n'], num2str(flmap,'%7.0f'));
fprintf(['ACCEPT m=4 only if its rms is comparable to m=2/m=3. Smaller-rms T4sgn\n' ...
         'is the correct sign. If m=2/m=3 themselves disagree badly, the HARNESS\n' ...
         'is still wrong -- do not read the m=4 column.\n']);
disp('MXFILL4_CHECK2_DONE');
