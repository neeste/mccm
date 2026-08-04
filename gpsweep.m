% GPSWEEP -- re-test gampro as a latency lever, at nch=3 with the margin s3fit built.
%
% WHY RE-OPEN IT. gampro was rejected 2026-08-01 as "stability-bankrupt": the
% exchange rate was ~2600 dB of oscillatory margin per unit of slope, and the
% nch=1 model had only 56.6 dB to spend, buying 0.021 against the 0.19 needed.
% THAT VERDICT WAS MEASURED AT nch=1 WITH 56.6 dB OF MARGIN. s3fit's nch=3 model
% sits at osc -180.2 dB -- 3.2x the margin -- and needs only 0.082 of slope
% (0.4950 -> 0.413). At the old rate that is ~213 dB, so it is close rather than
% hopeless, and the rate itself may differ at nch=3.
%
% latatlas also ranks gpgrade as the best real lever on the board: dslope
% +0.21428 with dlvlc only +0.6746 (S=0.32), against hbsc's S=0.48 at a
% twentieth the slope movement. It is the only parameter with BOTH large dslope
% and small level_c cost. Its cost is not level dependence, it is stability --
% which is exactly the column latatlas did not measure. This measures it.
%
% CF-MAP-FREE: maperr was EXACTLY constant across gampro grades at nch=1
% (90.97 at every g). Asserted here rather than assumed.
PROJ = '/Users/neely/Library/CloudStorage/OneDrive-FatherFlanagan''sBoysHome/STN/mod/mccm';
addpath(PROJ); cd(PROJ);
L=load('fit_nch3_surface.mat'); pa=L.R.pa;
if (isfield(pa,'hbmode')), pa=rmfield(pa,'hbmode'); end
H=parfit26('handles'); pv0=H.getpar(pa); ig=30+numel(pa.chsz)+1;
fprintf('\n  gpgrade base = %.6f  (pv index %d)\n', pv0(ig), ig);
gs=[0 -0.005 -0.010 -0.015 -0.020 -0.025 -0.030 -0.040];
fprintf('\n  %8s %8s %8s %9s %9s %9s\n','gpgrade','slope','lvl_c','shoulder','maperr','osc');
res=[];
for i=1:numel(gs)
    pvk=pv0; pvk(ig)=gs(i); p=H.setpar(pa,pvk,[]);
    try
        m=abr_metric(p,false);
        Rf=fdm26(struct('pa',p));
        evalc('C=tdm26(''coupeig'',struct(''pa'',p));');
        fprintf('  %8.4f %8.4f %8.3f %9.4f %9.2f %9.1f\n', ...
                gs(i), m.slope, m.level_c, m.shoulder, Rf.maperr, C.maxRe_osc);
        res(end+1,:)=[gs(i) m.slope m.level_c m.shoulder Rf.maperr C.maxRe_osc]; %#ok<SAGROW>
    catch e
        fprintf('  %8.4f  FAILED: %s\n', gs(i), e.message(1:min(50,end)));
    end
end
if (size(res,1)>1)
    fprintf('\n  maperr spread across the sweep: %.4f  (CF-map-free if ~0)\n', ...
            max(res(:,5))-min(res(:,5)));
    ok=res(res(:,6)<-40,:);
    if (~isempty(ok))
        [~,j]=min(abs(ok(:,2)-0.413));
        fprintf('  BEST SUB-CRITICAL: g %.4f -> slope %.4f (target 0.413), lvl_c %.3f, osc %.1f\n', ...
                ok(j,1), ok(j,2), ok(j,3), ok(j,6));
    end
end
save('/Users/neely/mccm_runs/gpsweep.mat','res','gs');
disp('GPSWEEP_DONE');
