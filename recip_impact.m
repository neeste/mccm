% How much does making the reciprocal stapes port the DEFAULT move the ABR
% metrics for the best previously-fitted parameter set?  Determines whether the
% prior fits need re-running under the corrected (bi-directional) physics.
L=load('parfit26_shoulder.mat'); R=L.R;
fn=fieldnames(R); fprintf('R fields: %s\n', strjoin(fn',', '));
pa=[];
if (isfield(R,'pa')), pa=R.pa; fprintf('using R.pa\n');
elseif (isfield(R,'pr')), pa=R.pr; fprintf('using R.pr\n');
elseif (isfield(R,'Rf') && isfield(R.Rf,'pa')), pa=R.Rf.pa; fprintf('using R.Rf.pa\n');
end
if (isempty(pa))
    fprintf('NO parameter struct stored in R -- cannot re-evaluate; listing nested:\n');
    for k=1:numel(fn)
        v=R.(fn{k});
        if (isstruct(v)), fprintf('  R.%s -> %s\n', fn{k}, strjoin(fieldnames(v)',', ')); end
    end
    return
end
fprintf('\nprev fit (legacy port): slope=%.3f level=%.3f anchor=%.2f shoulder=%.3f\n', ...
    R.mf.slope, 100*(R.mf.level_c^(1/100)-1), R.mf.lat(2,1), R.mf.shoulder);
for rc=[0 1]
    p=pa; p.me_recip=rc;
    m=abr_metric(p,false);
    if (~m.ok), fprintf('recip=%d : FAILED (%s)\n',rc,m.msg); continue; end
    fi=find(abs(m.f-1)<0.01,1); li=find(abs(m.slv-20)<0.1,1);
    fprintf('recip=%d : slope=%.3f  level=%.3f %%/dB  anchor=%.2f ms  shoulder=%.3f\n', ...
        rc, m.slope, 100*(m.level_c^(1/100)-1), m.lat(fi,li), m.shoulder);
end
disp('RECIP_IMPACT_DONE');
