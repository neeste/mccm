% Tuning/MAP comparison, 1- vs 3-chamber, alongside the corrected latency fit.
% maperr comes from fdm26 (frequency domain) -- unaffected by the tdm26 detector
% or the reciprocal-port change, so it is recomputed fresh and directly compared.
files={'parfit26_recip.mat','parfit26_shoulder.mat','parfit26_anchored.mat', ...
       'parfit26_corrected.mat','refit_c3_map.mat','refit_c1broad.mat'};
RB=[]; if (exist('rebaseline.mat','file')), L=load('rebaseline.mat'); RB=L.Res; end
fprintf('\n%-24s %4s | %9s %9s | %8s %7s\n','fit','nch','maperr','maperr','bandRMS','slope');
fprintf('%-24s %4s | %9s %9s | %8s %7s\n','','','(fresh)','(saved)','(corrected)','');
fprintf('%s\n',repmat('-',1,74));
for k=1:numel(files)
    fn=files{k};
    if (~exist(fn,'file')), fprintf('%-24s (missing)\n',fn); continue; end
    try
        L=load(fn); R=L.R; pa=R.pa;
        nch=1; if (isfield(pa,'m')), nch=pa.m; end
        sv=NaN; if (isfield(R,'Rf') && isfield(R.Rf,'maperr')), sv=R.Rf.maperr; end
        Rf=fdm26(struct('pa',pa));
        br=NaN; sl=NaN;
        for j=1:numel(RB)
            if (strcmp(RB(j).file,fn)), br=RB(j).rms; sl=RB(j).m.slope; end
        end
        fprintf('%-24s %4d | %9.1f %9.1f | %8.3f %7.3f\n', fn, nch, Rf.maperr, sv, br, sl);
    catch e
        fprintf('%-24s ERROR: %s\n', fn, e.message);
    end
end
fprintf('\n(lower maperr = better tuning/MAP fit; fits used maptol 167-185 as acceptance)\n');
disp('MAPERR_COMPARE_DONE');
