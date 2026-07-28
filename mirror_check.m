% VALIDATE THE fdm26 m3form MIRROR against tdm26.
%
% fdm26's m==3 zk now honours pa.m3form (it previously ignored it, so maperr was
% always computed with the d2-only law -- the reason nesting_test B and E gave
% identical maperr 499.3).
%
% THREE CHECKS, in order of what they can prove:
%  1 REGRESSION: with m3form=0 the maperr must be UNCHANGED (499.3 for m=3
%    native). If this moves, the edit broke the default path.
%  2 EFFECT: with m3form=1 the maperr must now DIFFER from 499.3. If it does
%    not, the switch still is not reaching the solve.
%  3 CORRECTNESS (the real one): fdm26 and tdm26 are independent implementations
%    of the same physics, so their CF maps must AGREE under BOTH settings. Peak
%    place per frequency, fdm26 via the cfmap dispatch vs tdm26 via the click --
%    the same cross-check that resolved T4sgn for the m=4 extension. m3form=0
%    CALIBRATES the comparison; only then is the m3form=1 agreement meaningful.

flmap = 500*2.^(-1:5);
ISV   = 1391:-10:11; XLO=0.05; XHI=0.85;

fprintf('\n=== CHECKS 1 & 2: maperr vs m3form (m=3 native) ===\n');
fprintf('  m3form | maperr   | verdict\n');
mp = nan(1,2);
for mf = [0 1]
    pa = modpar26(3); pa.m3form = mf;
    try
        R = fdm26(struct('pa',pa));
        mp(mf+1) = R.maperr;
        fprintf('    %d    | %8.1f |\n', mf, R.maperr);
    catch e
        fprintf('    %d    | FAILED: %s\n', mf, e.message);
    end
end
if (isfinite(mp(1)))
    if (abs(mp(1)-499.3) < 0.5), fprintf('  CHECK 1 PASS: default unchanged (499.3)\n');
    else, fprintf('  CHECK 1 *** FAIL ***: default moved 499.3 -> %.1f\n', mp(1)); end
end
if (all(isfinite(mp)))
    if (abs(mp(1)-mp(2)) > 0.5), fprintf('  CHECK 2 PASS: m3form now CHANGES maperr (delta %.1f)\n', abs(mp(1)-mp(2)));
    else, fprintf('  CHECK 2 *** FAIL ***: m3form still has no effect on maperr\n'); end
end

fprintf('\n=== CHECK 3: fdm26 vs tdm26 CF map (peak place x/L per frequency) ===\n');
for mf = [0 1]
    pa = modpar26(3); pa.m3form = mf;
    xf = nan(1,numel(flmap)); xt = xf;
    try
        Rc = fdm26(struct('cfmap',1,'pa',pa,'flst',flmap)); xf = Rc.xpk_bm;
    catch e
        fprintf('  m3form=%d cfmap FAILED: %s\n', mf, e.message); continue
    end
    try
        pt = pa; pt.isv = ISV;
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
        fprintf('  m3form=%d tdm26 FAILED: %s\n', mf, e.message); continue
    end
    d = xf - xt; rms = sqrt(mean(d(isfinite(d)).^2));
    fprintf('  m3form=%d  fdm26: %s\n', mf, num2str(xf,'%7.3f'));
    fprintf('            tdm26: %s   rms dx/L = %.4f\n', num2str(xt,'%7.3f'), rms);
end
fprintf('  freqs (Hz): %s\n', num2str(flmap,'%7.0f'));
fprintf(['\n  ACCEPT the mirror only if the m3form=1 rms is COMPARABLE to the\n' ...
         '  m3form=0 rms (which calibrates the comparison). A much larger rms at\n' ...
         '  m3form=1 means the frequency-domain form does not match tdm26.\n']);
disp('MIRROR_CHECK_DONE');
