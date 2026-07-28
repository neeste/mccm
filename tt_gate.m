% WHY does the tdm26 validity gate reject everything? Test tiptail_metric on the
% EXACT starting model of the fit (pinned m2o x32, native profiles) -- the check
% I should have run BEFORE launching (my own pre-flight rule).
b2=modpar26(4).m2o;
cfg={ {'m=2 native (gate reference)', modpar26(2)}, ...
      {'m=4 native',                  modpar26(4)}, ...
      {'m=4 m2o x32 (FIT START)',     setfield(modpar26(4),'m2o',b2*32)} };
for c=1:numel(cfg)
    nm=cfg{c}{1}; pa=cfg{c}{2};
    m=tiptail_metric(pa,false);
    fprintf('\n%s\n', nm);
    fprintf('   ok=%d  nvalid=%d  d=%.3f  R2=%.2f  msg=%s\n', m.ok, m.nvalid, m.d, m.r2, m.msg);
    if(~isempty(m.BF)), fprintf('   BF kept (kHz): %s\n', num2str(m.BF,'%7.2f')); end
end
fprintf('\nGATE: tiptail_metric needs >=4 places with BF in 0.6-16 kHz, strictly\n');
fprintf('increasing. It was tuned for the 2-chamber map (0.32-16 kHz). If the\n');
fprintf('4-chamber compressed map puts <4 places in that window, tt.ok=false and\n');
fprintf('jointobj returns J=1e6 for EVERY parameter set -> no descent possible.\n');
disp('TT_GATE_DONE');
