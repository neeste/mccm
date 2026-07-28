L=load('champions.mat'); C=L.C;
fprintf('banked configs: %d\n', numel(C));
for i=1:numel(C)
    S=C(i).S; fn=fieldnames(S);
    fprintf('\n m=%d  <- %s  (%.0f s)\n', C(i).nch, C(i).src, S.wall);
    if (isfield(S,'err')), fprintf('   ERROR: %s\n', S.err); continue; end
    fprintf('   TBABR block present: %d   TBOAE block present: %d\n', ...
            isfield(S,'surf_resid'), isfield(S,'oae_lat'));
    if (isfield(S,'tbabr_err')), fprintf('   TBABR FAILED: %s\n', S.tbabr_err); end
    if (isfield(S,'surf_err')),  fprintf('   SURF FAILED: %s\n', S.surf_err); end
    if (isfield(S,'surf_resid'))
        fprintf('   d=%.3f b=%.2f c=%.2f shape-rms=%.3f\n', S.d,S.b,S.c,S.surf_resid);
    end
    if (isfield(S,'oae_lat')), fprintf('   OAE lat=%.2f ms ratio=%.2f n=%d\n', S.oae_lat,S.oae_ratio,S.oae_n);
    elseif (isfield(S,'oae')), fprintf('   OAE present but all-NaN (n=%d finite)\n', nnz(isfinite(S.oae))); end
    fprintf('   maperr=%.1f maxRe=%+.1f range=%.2f fold=%.3f contrast=%.1f amp=%+.2f\n', ...
            S.maperr,S.maxRe,S.bf_range,S.bf_fold,S.contrast,S.amp_gain);
end
disp('PEEK_DONE');
