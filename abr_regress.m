function R = abr_regress(mode, varargin)
% ABR_REGRESS  Month-1 latency regression harness for the mccm cochlear model.
%   Scores tbabr WNR latency against the analytic ABR forward-latency target and
%   probes the basal cochlear-amplifier gain vs numerical-instability trade-off
%   that is shortening high-frequency latencies (MOH24 Fig 3 -> MATLAB).
%
%   R = ABR_REGRESS('baseline', NCH)
%       Run once for NCH chambers (default NCH=1). Prints a scorecard and plots
%       model-vs-target latency curves. NCH may be a vector, e.g. [1 3].
%
%   R = ABR_REGRESS('sweep', NCH, FIELD, VALUES)
%       Sweep a scalar model parameter FIELD (e.g. 'ace','aco','gam') over
%       VALUES. For each value: pa = modpar26(NCH); pa.(FIELD) = value.
%       Records error metrics and whether the run stayed stable, then plots
%       (a) error & slope vs FIELD with the stability boundary marked and
%       (b) latency-vs-frequency curves colored by FIELD value.
%       Example: abr_regress('sweep', 1, 'ace', -0.30:-0.05:-0.55)
%
%   R = ABR_REGRESS('gampro', NCH, BOOSTS, BASALFRAC)
%       Sweep a smooth basal boost of the named CA-gain profile pa.gampro.
%       gampro = 1 over most of the cochlea, ramped up by (1+BOOST) across the
%       basal fraction BASALFRAC (place from the stapes; default 0.35). BOOST=0
%       reproduces the baseline. This directly probes "restore basal gain".
%       Example: abr_regress('gampro', 1, [0 .5 1 2 4], 0.35)
%
%   All modes return R and save it to a .mat file. R.summary is a table-like
%   struct array (one row per condition) for quick inspection.
%
%   See also ABR_METRIC, TDM26, MODPAR26.

switch lower(mode)
    case 'baseline', R = do_baseline(varargin{:});
    case 'sweep',    R = do_sweep(varargin{:});
    case 'gampro',   R = do_gampro(varargin{:});
    otherwise, error('abr_regress: unknown mode "%s"', mode);
end
end

% ======================================================================
function R = do_baseline(nch)
if (nargin < 1), nch = 1; end
R.mode = 'baseline'; R.nch = nch;
% --- compute all metrics first (no figures open, so no figure(1) collision) ---
M = [];
for i = 1:numel(nch)
    fprintf('running baseline nch=%d (tbabr, ~%s) ...\n', nch(i), tern(nch(i)<3,'40s','95s'));
    m = abr_metric(nch(i), false);
    print_scorecard(m);
    M = [M; m]; %#ok<AGROW>
end
R.m = M;
% --- then plot ---
figure('Name','abr_regress baseline','Color','w'); clf
cols = lines(numel(nch));
leg = {};
for i = 1:numel(nch)
    m = M(i);
    if (~m.ok), continue; end
    subplot(1,2,1);
    loglog(m.f, m.lat, 'o-', 'Color', cols(i,:)); hold on
    loglog(m.f, m.abr, '--', 'Color', 0.6*[1 1 1]);
    subplot(1,2,2);
    loglog(m.f, m.rmse_f, 'o-', 'Color', cols(i,:)); hold on
    leg{end+1} = sprintf('nch=%d', m.nch); %#ok<AGROW>
end
subplot(1,2,1);
xlabel('frequency (kHz)'); ylabel('latency (ms)');
title('model (o-) vs target (--)'); axis([0.4 5 1 20]); grid on
subplot(1,2,2);
xlabel('frequency (kHz)'); ylabel('per-frequency RMSE (ms)');
title('where the error lives'); grid on; if(~isempty(leg)), legend(leg); end

fn = sprintf('abr_baseline_nch%s.mat', sprintf('_%d', nch));
save(fn, 'R'); fprintf('saved %s\n', fn);
end

% ======================================================================
function R = do_sweep(nch, field, values)
if (nargin < 3), error('usage: abr_regress(''sweep'', nch, field, values)'); end
values = values(:)';
R.mode = 'sweep'; R.nch = nch; R.field = field; R.values = values;
n = numel(values);
[rmse, logrmse, slope, ok, wall] = deal(nan(1,n));
lat = cell(1,n); msg = cell(1,n);
fprintf('\n== sweep nch=%d  %s over %d values ==\n', nch, field, n);
pa0 = get_base_pa(nch);
for j = 1:n
    pa = pa0;
    pa.(field) = values(j);
    m = abr_metric(pa, false);
    rmse(j)=m.rmse; logrmse(j)=m.logrmse; slope(j)=m.slope;
    ok(j)=m.ok; wall(j)=m.wall; lat{j}=m.lat; msg{j}=m.msg;
    R.m(j) = m;
    fprintf('%s=%+8.4g  %-6s  logRMSE=%6.3f  slope=%5.3f  c=%4.2f  (%4.1fs)  %s\n', ...
        field, values(j), tern(m.ok,'OK','UNSTBL'), logrmse(j), slope(j), m.level_c, wall(j), msg{j});
end
R.summary = struct('value',num2cell(values),'ok',num2cell(ok), ...
    'rmse',num2cell(rmse),'logrmse',num2cell(logrmse),'slope',num2cell(slope));
plot_sweep(R.m, field, values, logrmse, slope, ok, lat, nch);
fn = sprintf('abr_sweep_nch%d_%s.mat', nch, field);
save(fn, 'R'); fprintf('saved %s\n', fn);
end

% ======================================================================
function R = do_gampro(nch, boosts, basalfrac)
if (nargin < 3), basalfrac = 0.35; end
if (nargin < 2), boosts = [0 0.5 1 2 4]; end
boosts = boosts(:)';
R.mode = 'gampro'; R.nch = nch; R.boosts = boosts; R.basalfrac = basalfrac;
n = numel(boosts);
[logrmse, slope, ok, wall] = deal(nan(1,n));
lat = cell(1,n); msg = cell(1,n);
pa0 = get_base_pa(nch); N = pa0.n;
% smooth basal ramp: full boost at the stapes (index 1), tapering to 0 at basalfrac*L
xfrac = ((0:N-1)') / (N-1);            % 0 at base (stapes) .. 1 at apex
ramp  = 0.5*(1 + cos(pi * min(xfrac/basalfrac, 1)));  % 1 at base -> 0 beyond frac
fprintf('\n== gampro basal-boost sweep nch=%d  basalfrac=%.2f ==\n', nch, basalfrac);
for j = 1:n
    pa = pa0;
    pa.gampro = 1 + boosts(j) * ramp;   % length N, =1 apically, =1+boost basally
    m = abr_metric(pa, false);
    logrmse(j)=m.logrmse; slope(j)=m.slope; ok(j)=m.ok; wall(j)=m.wall;
    lat{j}=m.lat; msg{j}=m.msg; R.m(j) = m;
    fprintf('boost=%+6.2f  %-6s  logRMSE=%6.3f  slope=%5.3f  c=%4.2f  (%4.1fs)  %s\n', ...
        boosts(j), tern(m.ok,'OK','UNSTBL'), logrmse(j), slope(j), m.level_c, wall(j), msg{j});
end
R.summary = struct('boost',num2cell(boosts),'ok',num2cell(ok), ...
    'logrmse',num2cell(logrmse),'slope',num2cell(slope));
plot_sweep(R.m, 'gampro boost', boosts, logrmse, slope, ok, lat, nch);
fn = sprintf('abr_gampro_nch%d.mat', nch);
save(fn, 'R'); fprintf('saved %s\n', fn);
end

% ======================================================================
% helpers
% ======================================================================
function print_scorecard(m)
fprintf('\n--- abr_metric scorecard: nch=%d ---\n', m.nch);
if (~m.ok)
    fprintf('  UNSTABLE / invalid: %s  (%.1fs)\n', m.msg, m.wall); return
end
fprintf('  linear RMSE   : %6.3f ms\n', m.rmse);
fprintf('  log10  RMSE   : %6.3f  (octave-fair)\n', m.logrmse);
fprintf('  slope d       : %6.3f   (target %.3f; higher = too steep)\n', m.slope, m.slope_tgt);
fprintf('  level c       : %6.3f   (target %.1f; lower = level dep. too weak)\n', m.level_c, m.level_c_tgt);
fprintf('  per-freq RMSE : '); fprintf('%5.2f ', m.rmse_f); fprintf(' ms  @ %s kHz\n', mat2str(m.f'));
if (m.n_sub>0), fprintf('  sub-threshold : %d of %d conditions excluded (NaN latency)\n', m.n_sub, numel(m.lat)); end
fprintf('  wall time     : %.1f s\n', m.wall);
end

function plot_sweep(metrics, field, values, logrmse, slope, ok, lat, nch)
ok = logical(ok);
figure('Name',['abr_regress sweep: ' field],'Color','w'); clf
% recover the frequency grid and target from the first stable metric (no re-run)
f = []; abr = [];
for j = 1:numel(metrics)
    if (metrics(j).ok), f = metrics(j).f; abr = metrics(j).abr; break; end
end
% (a) error & slope vs swept value
subplot(1,2,1);
yyaxis left
plot(values, logrmse, 'o-'); ylabel('log10 RMSE (octave-fair)');
yyaxis right
plot(values, slope, 's--'); hold on
yline(0.413, ':', 'target slope');
ylabel('latency slope d');
xlabel(field);
title(sprintf('nch=%d: error & slope vs %s', nch, field));
grid on
% mark unstable points
bad = ~logical(ok);
if (any(bad) && any(ok))
    yyaxis left; hold on
    ybad = min(logrmse(ok),[],'omitnan');
    plot(values(bad), zeros(1,nnz(bad)) + ybad, 'rx', 'MarkerSize', 10, 'LineWidth', 1.5);
    text(values(find(bad,1)), double(ybad), '  unstable', 'Color','r');
end
% (b) latency-vs-frequency family, colored by swept value
subplot(1,2,2);
cols = parula(numel(values));
leg = {};
for j = 1:numel(values)
    if (~ok(j) || isempty(lat{j})), continue; end
    loglog(f, mean(lat{j},2), 'o-', 'Color', cols(j,:)); hold on
    leg{end+1} = sprintf('%s=%.3g', field, values(j)); %#ok<AGROW>
end
if (~isempty(abr))
    loglog(f, mean(abr,2), 'k--', 'LineWidth', 1.5);
    leg{end+1} = 'target';
end
xlabel('frequency (kHz)'); ylabel('level-mean latency (ms)');
axis([0.4 5 1 20]); grid on;
title('latency family (level-averaged)');
if (~isempty(leg)), legend(leg, 'Location','southwest'); end
end

function s = tern(c, a, b), if c, s = a; else, s = b; end, end

function pa = get_base_pa(nch)
% Return a base parameter struct for NCH chambers. modpar26 is now a shared
% standalone file (modpar26.m), so call it directly.
pa = modpar26(nch);
end
