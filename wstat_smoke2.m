% Does wstat actually CHANGE J? The first smoke test did not test this: m=1's
% maxRe sat below statol in both arms, so the static term was zero either way
% and the runs were identical for a reason unrelated to wstat.
% Here statol is driven far negative so max(0, maxRe - statol) is guaranteed
% LARGE and active, and only wstat differs. The START J (iter 0) must scale.
pa = modpar26(1);
b = struct('wslope',0,'wlevel',0,'wlcb',0,'skipabr',1,'cheapstab',0, ...
           'maptol',0,'fitidx',[1 2 3],'warm',pa,'maxfe',1, ...
           'pin',struct('chsz',[1 1]),'statol',-1000,'out','wstat_smoke2.mat');
for w = [0 0.001 0.01]
    o = b; o.wstat = w;
    fprintf('\n  ===== wstat = %g (statol %g, term forced active) =====\n', w, o.statol);
    parfit26(1, o);
end
fprintf(['\n  The iter-0 J values above must DIFFER and scale with wstat. If they\n' ...
         '  are identical the option is inert and the edit did not take.\n']);
disp('WSTAT_SMOKE2_DONE');
