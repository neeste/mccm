% What does the click/impulse path return, and does it run for m=2,3,4?
% modpar26 sets hbnl=0 by default, so the click response is already LINEAR
% (tbabr_condition is what sets hbnl=1). stcfg=0 -> short_chirp (broadband).
for nch=[2 3 4]
    fprintf('\n--- nch=%d ---\n',nch);
    try
        pa=modpar26(nch);
        fprintf('  params: parlab=%s  m=%d  hbnl=%d  chsz=[%s]  isv n=%d\n', ...
                pa.parlab, pa.m, pa.hbnl, num2str(pa.chsz), numel(pa.isv));
        t0=tic; S=tdm26(0,nch,0,0); w=toc(t0);
        fn=fieldnames(S);
        fprintf('  ran OK in %.0f s;  returned %d fields\n', w, numel(fn));
        has=@(f) any(strcmp(fn,f));
        fprintf('  has bf=%d mbf=%d d1=%d d2=%d f=%d\n', has('bf'),has('mbf'),has('d1'),has('d2'),has('f'));
        if (has('bf')),  fprintf('  bf  = %s\n', num2str(S.bf,'%8.3f')); end
        if (has('mbf')), fprintf('  mbf = %s\n', num2str(S.mbf,'%8.3f')); end
        if (has('d1')),  fprintf('  d1 size = %s   f range = %.3f..%.1f kHz\n', mat2str(size(S.d1)), min(S.f), max(S.f)); end
    catch e
        fprintf('  FAILED: %s\n', e.message);
        if (~isempty(e.stack)), fprintf('    at %s line %d\n', e.stack(1).name, e.stack(1).line); end
    end
end
disp('TUNE_PROBE_DONE');
