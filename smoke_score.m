% Pre-flight: does score26 run, and is the fast block sane on a known model?
% m=2 native is the control -- its map should span ~5.9 oct, monotonic.
S=score26(modpar26(2),'fast');
fprintf('\nCONTROL CHECK (m=2 native): range %.2f oct (expect ~5.9), mono %s, maperr %.1f\n', ...
        S.bf_range, S.bf_mono, S.maperr);
disp('SMOKE_SCORE_DONE');
