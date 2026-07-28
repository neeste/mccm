% VALIDATE THE m=4 EXTENSION OF mxfill AGAINST tdm26.
%
% fdm26 (frequency domain) and tdm26 (time domain) are INDEPENDENT
% implementations of the same physics, so their CF maps must agree. Strategy:
%   1. calibrate the comparison on m=2 and m=3, where fdm26 is long-established;
%   2. apply the SAME comparison to m=4 and require comparable agreement;
%   3. resolve the one free sign in T (pa.T4sgn = +1 / -1) empirically -- the
%      derivation fixes T up to that sign, so tdm26 decides it, not a guess.
%
% Observable: BF per saved place. tdm26 gives it from the linear click; fdm26
% gives it from the peak of |Ybm| vs frequency at the same place. Both are the
% place-frequency map, which is the thing any tuning/maperr number rests on.

fr = 2.^(linspace(log2(0.5), log2(12), 40));   % kHz probe grid for fdm26
ISV = [1136 1005 840 655 466 273 80];          % the standard saved places

fprintf('\n=== fdm26 runs at all? ===\n');
for nch=[2 3 4]
    for sg=[1 -1]
        pa=modpar26(nch); if (nch==4), pa.T4sgn=sg; end
        try
            R=fdm26(struct('pa',pa));
            me='n/a'; if (isstruct(R)&&isfield(R,'maperr')), me=sprintf('%.2f',R.maperr); end
            fprintf('  nch=%d T4sgn=%+d : OK   maperr=%s\n', nch, sg, me);
        catch e
            fprintf('  nch=%d T4sgn=%+d : THREW -> %s\n', nch, sg, e.message);
        end
        if (nch~=4), break; end
    end
end

fprintf('\n=== CF-map agreement: tdm26 (click) vs fdm26 (|Ybm| peak) ===\n');
fprintf('  nch  T4sgn |      BF per place (kHz), tdm26 over fdm26      | rms log2 diff\n');
for nch=[2 3 4]
  sgl = 1; if (nch==4), sgl=[1 -1]; end
  for sg=sgl
    pa=modpar26(nch); pa.isv=ISV; if (nch==4), pa.T4sgn=sg; end
    bt=nan(1,numel(ISV)); bf=bt;
    % --- tdm26 CF map from the linear click ---
    try
        evalc('S=tdm26(0,pa,0,0);');
        nf=numel(S.f); f=S.f(:); P=fft(S.ped); P=P(1:nf);
        for i=1:numel(ISV)
            D=fft(S.d1(:,i)); H=D(1:nf)./max(abs(P),eps);
            a=abs(H).*f; a(~(f>0.3&f<16))=0;
            [pk,ip]=max(a); if (pk>0), bt(i)=f(ip); end
        end
    catch e
        fprintf('  nch=%d: tdm26 failed: %s\n', nch, e.message); continue
    end
    % --- fdm26 CF map: peak of |Ybm| vs frequency at each place ---
    try
        amp=zeros(numel(fr),numel(ISV));
        for j=1:numel(fr)
            Rj=fdm26(struct('fdsolve',1,'pa',pa,'f',fr(j)*1000));
            yb=[];
            if (isstruct(Rj) && isfield(Rj,'Yb')), yb=Rj.Yb; end   % fd_solve returns Yb
            if (isempty(yb)), error('no Yb field from fdsolve'); end
            yb=yb(:,1);                                  % BM column if multi-column
            amp(j,:)=abs(yb(ISV));
        end
        for i=1:numel(ISV)
            [pk,ip]=max(amp(:,i)); if (pk>0), bf(i)=fr(ip); end
        end
    catch e
        fprintf('  nch=%d T4sgn=%+d: fdm26 CF map failed: %s\n', nch, sg, e.message);
        fprintf('        tdm26 BF: %s\n', num2str(bt,'%7.2f')); continue
    end
    d=log2(bt./bf); rms=sqrt(mean(d(isfinite(d)).^2));
    fprintf('  %d   %+5d | t: %s\n', nch, sg, num2str(bt,'%7.2f'));
    fprintf('              | f: %s |  %.3f oct\n', num2str(bf,'%7.2f'), rms);
  end
end
fprintf(['\nACCEPT the m=4 extension only if its rms log2 CF-map difference is\n' ...
         'COMPARABLE to m=2/m=3 (which validate the comparison itself). The T4sgn\n' ...
         'with the smaller rms is the correct sign. A large rms for BOTH signs\n' ...
         'means the extension is wrong somewhere else (zk, D, or B), not just T.\n']);
disp('MXFILL4_CHECK_DONE');
