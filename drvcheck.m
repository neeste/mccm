% DRVCHECK -- is the c3fit slope an artifact of the latency drive?
%
% parfit26.m:183 carries a standing warning: "latency drive ('bm' = BM-peak).
% SET THIS for valid fits: the default drive rewards the WNR peak-switch
% artifact (that is what produced the bogus 0.413)." It predates the git
% history, so it cannot be dated -- only measured.
%
% NEITHER longfit NOR c3fit SET opts.hbmode, and neither set wshoulder, so both
% guards were off for every slope on record. If the slope moves materially
% under hbmode='bm', the 0.413 results are drive-dependent and s3fit must pin
% the drive. Cheap to settle now that abr_metric is ~7 s.
%
% cur.hb = v1 - v2 (shear) by default; 'bm' = v1 alone (tdm26 hc_step).
PROJ = '/Users/neely/Library/CloudStorage/OneDrive-FatherFlanagan''sBoysHome/STN/mod/mccm';
addpath(PROJ);
S3 = load('/Users/neely/mccm_runs/c3fit_seg3.mat'); pa = S3.R.pa;

cfg = {'', 'bm', 'dof2'};
R = struct('mode',{},'slope',{},'lvlc',{},'shoulder',{},'nsub',{},'rmse',{});
M = cell(1,3);
for i = 1:numel(cfg)
    p = pa;
    if (isempty(cfg{i})), p = rmfield_safe(p,'hbmode'); else, p.hbmode = cfg{i}; end
    t=tic; m = abr_metric(p,false); el=toc(t);
    M{i} = m;
    R(end+1) = struct('mode',ternlab(cfg{i}),'slope',m.slope,'lvlc',m.level_c, ...
        'shoulder',m.shoulder,'nsub',m.n_sub,'rmse',m.rmse); %#ok<SAGROW>
    fprintf('  %-8s slope %.4f | lvl_c %6.3f | shoulder %.4f | n_sub %d | rmse %.3f | %.1f s\n', ...
            ternlab(cfg{i}), m.slope, m.level_c, m.shoulder, m.n_sub, m.rmse, el);
end

fprintf('\n=== DRIVE DEPENDENCE OF THE c3fit RESULT ===\n');
fprintf('  target slope 0.413; c3fit seg3 reported 0.4118 under the DEFAULT drive\n\n');
fprintf('  %-8s %8s %8s %10s\n','drive','slope','lvl_c','shoulder');
for i=1:numel(R)
    fprintf('  %-8s %8.4f %8.3f %10.4f\n', R(i).mode, R(i).slope, R(i).lvlc, R(i).shoulder);
end
fprintf('\n  shoulder: 0 = clean single-peaked CAP, ->1 = a competing 2nd peak\n');
fprintf('  |slope(bm) - slope(default)| = %.4f  (gap to target from default = %.4f)\n', ...
        abs(R(2).slope-R(1).slope), abs(R(1).slope-0.413));

fprintf('\n=== per-cell shoulder, DEFAULT drive (rows = f 0.5/1/2/4 kHz, cols = 20/40/60/80 dB) ===\n');
disp(round(M{1}.sho,3));
fprintf('=== per-cell latency (ms), DEFAULT drive ===\n');
disp(round(M{1}.lat,2));
fprintf('=== per-cell latency (ms), hbmode=bm ===\n');
disp(round(M{2}.lat,2));

save('/Users/neely/mccm_runs/drvcheck.mat','R','M','pa');
disp('DRVCHECK_DONE');

function p=rmfield_safe(p,f), if (isfield(p,f)), p=rmfield(p,f); end, end
function s=ternlab(c), if isempty(c), s='default'; else, s=c; end, end
