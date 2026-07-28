m=checkcode('tdm26.m','-struct');
syn=m(arrayfun(@(z)any(cellfun(@(s)~isempty(strfind(lower(z.message),s)),{'parse','syntax','unbalanced'})),m));
if (~isempty(syn))
    for k=1:numel(syn), fprintf('SYNTAX L%d %s\n',syn(k).line,syn(k).message); end
    return
end
fprintf('checkcode OK\n');
pr.fr=2; pr.lv=60; pr.pa=modpar26(3); pr.pa.oae=1; pr.pa.rough_amp=3e-2;
Sd=tdm26('wnr1',pr,0,0);                      % DEFAULT (no me_recip field)
pr.pa.me_recip=1; S1=tdm26('wnr1',pr,0,0);    % explicit reciprocal
pr.pa.me_recip=0; S0=tdm26('wnr1',pr,0,0);    % explicit legacy
fprintf('2kHz/60  default : WNR=%.4f  emag=%.1f\n', Sd.tpk, Sd.oam);
fprintf('         recip=1 : WNR=%.4f  emag=%.1f   (must MATCH default)\n', S1.tpk, S1.oam);
fprintf('         recip=0 : WNR=%.4f  emag=%.1f   (legacy, must differ)\n', S0.tpk, S0.oam);
fprintf('default==recip1 : %d\n', abs(Sd.tpk-S1.tpk)<1e-9 && abs(Sd.oam-S1.oam)<1e-9);
fprintf('default~=recip0 : %d\n', abs(Sd.tpk-S0.tpk)>1e-9 || abs(Sd.oam-S0.oam)>1e-9);
disp('VERIFY_DEFAULT_DONE');
