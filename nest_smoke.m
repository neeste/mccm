% Does the nested topology run, and does the old path still reproduce exactly?
b2=modpar26(4).m2o;
base=modpar26(4); base.m2o=b2*32;
cfg={ {'current (nested=0)', 0, 0}, {'nested, sealed (clvent=0)', 1, 0}, ...
      {'nested, vent 1', 1, 1}, {'nested, vent 10', 1, 10} };
fprintf('\n  config                     | d1 nonfin | max|d1|    | ped max   | BF apex..base (kHz)\n');
fprintf('%s\n',repmat('-',1,92));
for c=1:numel(cfg)
    nm=cfg{c}{1}; pa=base; pa.nested=cfg{c}{2}; pa.clvent=cfg{c}{3};
    pa.isv=[1136 1005 840 655 466 273 80];
    err='';
    try, evalc('S=tdm26(0,pa,0,0);'); catch e, err=e.message; end
    if(~isempty(err)), fprintf('  %-26s | THREW: %s\n', nm, err); continue; end
    nf=sum(~isfinite(S.d1(:)));
    nfd=numel(S.f); f=S.f(:); P=fft(S.ped); P=P(1:nfd); bf=nan(1,7);
    for i=1:7
        D=fft(S.d1(:,i)); H=D(1:nfd)./max(abs(P),eps); a=abs(H).*f; a(~(f>0.15&f<18))=0;
        [pk,ip]=max(a); if(pk>0), bf(i)=f(ip); end
    end
    fprintf('  %-26s | %9d | %10.3e | %9.3e | %s\n', nm, nf, ...
        max(abs(S.d1(isfinite(S.d1)))), max(abs(S.ped)), num2str(bf,'%6.2f'));
end
fprintf('\n  Row 1 must match the pre-change 4-chamber exactly (regression).\n');
fprintf('  Rows 2-4 are the nested chain; the vent should RAISE max|d1| toward\n');
fprintf('  the 3-chamber scale if the sealed-CL loading was the problem.\n');
disp('NEST_SMOKE_DONE');
