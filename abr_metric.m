function m = abr_metric(pr, verbose)
% ABR_METRIC  Run the toneburst-ABR (tbabr) protocol and score the model's
%   whole-nerve-response (WNR) latency against the analytic ABR forward-latency
%   target. This is the core, reusable metric for the Month-1 latency regression.
%
%   m = ABR_METRIC(NCH)         run default parameters for NCH fluid chambers.
%   m = ABR_METRIC(PA)          run with a parameter struct (from MODPAR26 with
%                               any fields, e.g. aco/ace/gampro, modified).
%   m = ABR_METRIC(..., false)  suppress tbabr console output (default true).
%
%   The metric wraps tdm26('tbabr',PR), which returns model latency S.lat and
%   the analytic target S.abr = b*c^(-lv/100)*fr^(-d) on a 4-freq x 4-level grid
%   (f = 0.5,1,2,4 kHz; lv = 20,40,60,80 dB SPL).
%
%   Returns struct M:
%     .ok        true if latencies are finite and physically plausible (a
%                divergence / numerical-instability check for gain sweeps)
%     .msg       reason string when .ok is false
%     .f .slv    frequency (kHz) and level (dB SPL) grids
%     .lat .abr  model and target latency (ms), 4x4 (freq x level)
%     .lev       achieved eardrum level (dB SPL), 4x4
%     .rmse      RMS(lat-abr) over all 16 conditions, ms   [linear, legacy]
%     .logrmse   RMS(log10 lat - log10 abr)                [octave-fair]
%     .rmse_f    per-frequency RMSE (4x1); rmse_f(4) is the basal (4 kHz) error
%     .slope     model latency slope d in  lat ~ f^(-d), averaged over level
%     .slope_tgt target slope (0.413); slope > slope_tgt means "too steep"
%     .level_c   model level dependence as effective c in lat ~ c^(-lv/100),
%                averaged over frequency; target 5. level_c < 5 => level
%                dependence too weak (latency doesn't shorten enough with level)
%     .level_c_f per-frequency effective c (4x1)
%     .level_c_tgt target c (5)
%     .nch .wall number of chambers, wall-clock seconds
%
%   See also ABR_REGRESS, TDM26, MODPAR26.

if (nargin < 2), verbose = true; end

if (isstruct(pr)), nch = pr.m; else, nch = max(1, pr); end

m.nch = nch;
m.slope_tgt = 0.413;   % d in the analytic target lat = b*c^(-i)*fr^(-d)
m.ok = false; m.msg = '';
m.rmse = NaN; m.logrmse = NaN; m.slope = NaN;
m.level_c = NaN; m.level_c_f = nan(4,1); m.level_c_tgt = 5;
m.n_sub = 0;   % count of sub-threshold (NaN-latency) conditions excluded
m.rmse_f = nan(4,1); m.lat = []; m.abr = []; m.lev = []; m.f = []; m.slv = [];
m.shoulder = NaN; m.sho = [];   % WNR shoulder ratio (2nd peak/main), mean over supra-threshold
m.oae = []; m.oam = [];         % tone-burst OAE latency/level grids (populated only if pr.oae set)

% Pass dsp1=dsp2=0 (tdm26 args 3&4) so the model suppresses ALL figure
% generation -- including its per-step animation (plot_ts -> drawnow), which
% otherwise renders ~1e4 times and turns a ~40 s run into hours interactively.
figs0 = findobj('Type','figure');   % belt-and-suspenders: close any stray figure
t0 = tic;
try
    if (verbose)
        S = tdm26('tbabr', pr, 0, 0);
    else
        evalc('S = tdm26(''tbabr'', pr, 0, 0);');   % run quietly, no figures
    end
catch e
    m.wall = toc(t0);
    m.msg = sprintf('tdm26 threw: %s', e.message);
    close_new_figs(figs0);
    return
end
m.wall = toc(t0);
close_new_figs(figs0);            % discard the auto tbabr figure; we re-plot

% ---- unpack ----
m.f = S.f(:); m.slv = S.slv(:)'; m.lat = S.lat; m.abr = S.abr; m.lev = S.lev;
if (isfield(S,'sho')), m.sho = S.sho; end
if (isfield(S,'oae')), m.oae = S.oae; end
if (isfield(S,'oam')), m.oam = S.oam; end

% ---- stability / plausibility check ----
% NaN latency = sub-threshold condition (the detector found no onset); exclude it.
% Inf or finite-but-implausible latency, or a bad achieved level = divergence: fail.
lat = m.lat;
fin = isfinite(lat);
m.n_sub = nnz(~fin);
if (any(~isfinite(m.lev(:))) || any(abs(m.lev - m.slv) > 25, 'all'))
    m.msg = 'achieved level far from requested (unstable middle-ear/cochlea)'; return
elseif (any(lat(fin) <= 0.2) || any(lat(fin) >= 40))
    m.msg = 'finite latency outside plausible [0.2, 40] ms band (diverged)'; return
elseif (nnz(fin) < numel(lat)/2)
    m.msg = sprintf('%d/%d conditions sub-threshold; too few to score', m.n_sub, numel(lat)); return
end
m.ok = true;
if (~isempty(m.sho) && any(fin(:))), m.shoulder = mean(m.sho(fin)); end   % mean shoulder, supra-threshold

% ---- error metrics over supra-threshold (finite) conditions only ----
d = lat - m.abr;
m.rmse    = sqrt(mean(d(fin).^2));
ld = log10(lat) - log10(m.abr);
m.logrmse = sqrt(mean(ld(fin).^2));
m.rmse_f  = nan(numel(m.f),1);
for j = 1:numel(m.f)
    fj = fin(j,:); if (any(fj)), m.rmse_f(j) = sqrt(mean(d(j,fj).^2)); end
end

% ---- slope of model latency vs frequency (per level, over finite freqs) ----
lf = log10(m.f);
dd = [];
for k = 1:numel(m.slv)
    fk = fin(:,k);
    if (nnz(fk) >= 2)
        p = polyfit(lf(fk), log10(lat(fk,k)), 1);
        dd(end+1) = -p(1); %#ok<AGROW>   lat ~ f^(-d) => d = -slope
    end
end
if (~isempty(dd)), m.slope = mean(dd); end

% ---- level slope: effective c (per freq, over finite levels) ----
% target: lat = b*c^(-lv/100), so log10(lat) is linear in level with slope
% -log10(c)/100. Fit per frequency and convert to an effective c (target 5).
lv = m.slv(:);
cf = [];
for j = 1:numel(m.f)
    fj = fin(j,:)';
    if (nnz(fj) >= 2)
        p = polyfit(lv(fj), log10(lat(j,fj))', 1);
        m.level_c_f(j) = 10^(-100 * p(1));
        cf(end+1) = m.level_c_f(j); %#ok<AGROW>
    end
end
if (~isempty(cf)), m.level_c = mean(cf); end
end

function close_new_figs(figs0)
newf = setdiff(findobj('Type','figure'), figs0);
if (~isempty(newf)), close(newf); end
end
