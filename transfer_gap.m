% THE LINEARIZED-vs-TONE-BURST TRANSFER GAP.
%
% Symptom: the 2-chamber tiptail fit drove the LINEARIZED (click) exponent to
% d=0.400 -- exactly on target -- while the TONE-BURST ABR slope stayed at
% 0.686. Same model, two measurements, 0.29 apart.
%
% HYPOTHESIS (from SN's domain point that the frequency-place map is LEVEL
% DEPENDENT -- the excitation maximum shifts with level): the click is a
% SMALL-SIGNAL measurement, i.e. the LOW-LEVEL limit, whereas the tone-burst
% slope reported by abr_metric is the MEAN of per-level slopes over 20-80 dB
% (abr_metric.m:104-109 fits log10(lat) vs log10(f) separately at each level and
% averages). If d rises with level, the mean is pulled far above the low-level
% value and the two measurements CANNOT agree -- not a bug, but two different
% quantities.
%
% TEST 1 is FREE: champions.mat already holds the full 4x4 latency surface for
% all four chamber counts, so the per-level slopes need no new simulation.
%   PREDICTION: d(20 dB) ~ the linearized value, rising with level.
%   REFUTED IF: d is flat across level -- then the gap is NOT level dependence
%   and the cause is the observable (BM group delay vs WNR peak) or the
%   nonlinearity (click runs hbnl=0, tbabr runs hbnl=1).
% TEST 2 runs tiptail_metric on the same champions for the linearized d.

L = load('champions.mat'); C = L.C;

fprintf('\n=== TEST 1: per-LEVEL tone-burst slope d (free, from champions.mat) ===\n');
fprintf('  m |');
fprintf('  d@%2.0fdB |', [20 40 60 80]); fprintf('  mean  | spread\n');
D = nan(numel(C),4);
for i = 1:numel(C)
    S = C(i).S;
    if (~isfield(S,'lat') || isempty(S.lat)), fprintf('  %d | (no latency surface)\n', C(i).nch); continue; end
    lat = S.lat; f = S.f(:); slv = S.slv(:);
    lf = log10(f);
    for k = 1:numel(slv)
        col = lat(:,k); ok = isfinite(col) & col > 0;
        if (nnz(ok) >= 2)
            p = polyfit(lf(ok), log10(col(ok)), 1);
            D(i,k) = -p(1);
        end
    end
    fprintf('  %d |', C(i).nch);
    fprintf('  %6.3f |', D(i,:));
    fprintf('  %5.3f | %5.3f\n', mean(D(i,isfinite(D(i,:)))), ...
            max(D(i,:))-min(D(i,:)));
end
fprintf('  levels are %s dB\n', num2str(C(1).S.slv));

fprintf('\n=== TEST 2: linearized (click) d on the SAME champions ===\n');
fprintf('  m | linearized d | R2   | nvalid | note\n');
for i = 1:numel(C)
    S = C(i).S;
    pa = get_pa(C(i).nch);
    if (isempty(pa)), fprintf('  %d | (params unavailable)\n', C(i).nch); continue; end
    try
        t = tiptail_metric(pa, false);
        if (t.ok || t.nvalid >= 4)
            fprintf('  %d |   %8.3f   | %4.2f |   %d    | %s\n', C(i).nch, t.d, t.r2, t.nvalid, t.msg);
        else
            fprintf('  %d |      --      |  --  |   %d    | %s\n', C(i).nch, t.nvalid, t.msg);
        end
    catch e
        fprintf('  %d | FAILED: %s\n', C(i).nch, e.message);
    end
end

fprintf(['\nREAD: if d(20 dB) is close to the linearized d and rises steeply with\n' ...
         'level, the transfer gap IS the level dependence -- the click measures the\n' ...
         'low-level limit and the tone-burst reports a level-average, so fitting one\n' ...
         'cannot move the other. If d is FLAT across level, the level hypothesis is\n' ...
         'refuted and the cause is the observable or the nonlinearity.\n']);
disp('TRANSFER_GAP_DONE');

function pa = get_pa(nch)
pa = [];
switch nch
    case 1
        if (exist('refit_c1broad.mat','file')), L=load('refit_c1broad.mat'); pa=L.R.pa; else, pa=modpar26(1); end
    case 2
        if (exist('parfit26_lowdamp.mat','file')), L=load('parfit26_lowdamp.mat'); pa=L.R.pa; else, pa=modpar26(2); end
    case 3
        if (exist('refit_c3_map.mat','file')), L=load('refit_c3_map.mat'); pa=L.R.pa; else, pa=modpar26(3); end
    case 4
        pa = modpar26(4); pa.m2o = pa.m2o*32;
end
end
