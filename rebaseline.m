% Re-baseline every milestone fit under the CORRECTED metric:
%   - reciprocal middle-ear port (now the default), and
%   - BM-peak latency measurement (hbmode='bm', verified non-perturbing).
% Reports the surface fit (band-RMS vs the 1988/2013 band, free Delta) and the
% worst per-20 dB step, which quantifies the compression cliff.
files={'parfit26_recip.mat','parfit26_shoulder.mat','parfit26_anchored.mat', ...
       'parfit26_corrected.mat','refit_c3_map.mat','refit_c1broad.mat'};
fprintf('\n%-24s | %-18s | %-34s | %s\n','fit','as-saved (old)','BM-peak metric (corrected)','steps');
fprintf('%-24s | %-18s | %-34s | %s\n','','slope  level_c','slope  level_c  bandRMS  Delta','worst per-20dB');
fprintf('%s\n',repmat('-',1,110));
Res=struct([]);
for k=1:numel(files)
    fn=files{k};
    if (~exist(fn,'file')), fprintf('%-24s | (missing)\n',fn); continue; end
    try
        L=load(fn); R=L.R; pa=R.pa; pa.hbmode='bm';
        os=NaN; ol=NaN;
        if (isfield(R,'mf')), if(isfield(R.mf,'slope')),os=R.mf.slope; end; if(isfield(R.mf,'level_c')),ol=R.mf.level_c; end; end
        m=abr_metric(pa,false);
        if (~m.ok), fprintf('%-24s | %5.3f %6.2f    | FAILED: %s\n',fn,os,ol,m.msg); continue; end
        f=m.f(:); slv=m.slv(:); [F,I]=ndgrid(f,slv/100);
        t88=12.90*(5.00.^(-I)).*(F.^(-0.413)); t13=12.63*(5.34.^(-I)).*(F.^(-0.390));
        tlo=min(t88,t13); thi=max(t88,t13);
        Lm=m.lat; ok=isfinite(Lm)&Lm>0;
        ds=linspace(0,1.5,301); best=inf; D=0;
        for j=1:numel(ds)
            V=max(0,tlo-(Lm+ds(j)))+max(0,(Lm+ds(j))-thi);
            r=sqrt(mean(V(ok).^2)); if(r<best),best=r; D=ds(j); end
        end
        st=0;
        for a=1:numel(f)
            for b=1:numel(slv)-1
                if (ok(a,b)&&ok(a,b+1))
                    st=max(st, abs(100*(Lm(a,b+1)-Lm(a,b))/Lm(a,b)));
                end
            end
        end
        fprintf('%-24s | %5.3f %6.2f    | %5.3f %6.2f  %7.3f  %5.2f    | %5.1f%%\n', ...
                fn, os, ol, m.slope, m.level_c, best, D, st);
        Res(end+1).file=fn; Res(end).m=m; Res(end).rms=best; Res(end).delta=D; Res(end).step=st; %#ok<SAGROW>
    catch e
        fprintf('%-24s | ERROR: %s\n', fn, e.message);
    end
end
fprintf('\n(target: slope 0.39-0.41, level_c 5.0-5.34, bandRMS ->0, steps ~28%%/20dB uniform)\n');
save('rebaseline.mat','Res');
disp('REBASELINE_DONE');
