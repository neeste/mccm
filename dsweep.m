% DSWEEP -- two diagnostics on the finished d3fit arms.
%
% 1. THE NEURAL-DELAY BOUND. The surface objective scores the 16-cell latency
%    matrix against a band, with Delta (neural delay in EXCESS of the assumed
%    5 ms) free in [0, 1.5] ms. Arm A returned Delta = 1.50 in BOTH segments,
%    i.e. AT the bound, while arm B returned 0.60 with room to spare. A term
%    reported under a binding constraint understates the disagreement, so arm
%    A's 1.342 ms is not comparable to arm B's 0.536 ms at face value. Sweeping
%    the ceiling says how much of arm A's misfit the bound is hiding, and
%    whether it is hiding any of arm B's.
%
%    Reads the SHIPPED surf_term through the handles hook rather than
%    reimplementing the band loss. A second copy of the objective is exactly the
%    duplication this project's standing rule forbids, and it is how the wlcb
%    band term went missing from a reported result once already.
%
% 2. THE SHOULDER, CELL BY CELL. The scalar shoulder is a mean over 16 cells and
%    can hide its own distribution: a mean of 0.26 is one thing if it is 16 cells
%    near 0.26 and quite another if it is 12 cells at 0.0 and 4 at 1.0, because
%    the detector only switches where the two peaks are COMPARABLE. Arm A's mean
%    is 0.0022 and arm B's is 0.2603; this prints the matrices so the difference
%    can be read as a distribution rather than as a summary statistic.
PROJ = '/Users/neely/Library/CloudStorage/OneDrive-FatherFlanagan''sBoysHome/STN/mod/mccm';
addpath(PROJ); cd(PROJ);
RUNDIR = '/Users/neely/mccm_runs';
H = parfit26('handles');

% arm B seg 2 may or may not have landed; prefer it, fall back to seg 1
fB = fullfile(RUNDIR,'d3fit_B_seg2.mat');
if (~exist(fB,'file')), fB = fullfile(RUNDIR,'d3fit_B_seg1.mat'); end
arms = {'A (d3int=1)', fullfile(RUNDIR,'d3fit_A_seg2.mat'); ...
        'B (control)', fB};

M = {};
for a=1:size(arms,1)
    L=load(arms{a,2}); pa=L.R.pa;
    if (isfield(pa,'hbmode')), pa=rmfield(pa,'hbmode'); end
    M{a} = abr_metric(pa,false); %#ok<SAGROW>
    fprintf('\n  %s   <- %s\n', arms{a,1}, arms{a,2});
end

% ---- 1. Delta-ceiling sweep -------------------------------------------------
fprintf('\n===== NEURAL-DELAY BOUND: band-RMS (ms) vs the ceiling on Delta =====\n');
fprintf('  %8s | %-24s | %-24s\n', 'ceiling', arms{1,1}, arms{2,1});
fprintf('  %8s | %10s %12s | %10s %12s\n', '(ms)', 'band-RMS', 'Delta used', 'band-RMS', 'Delta used');
for dhi = [1.5 2.0 2.5 3.0 4.0 5.0 6.0 8.0]
    P.dlo = 0; P.dhi = dhi;
    [s1,d1] = H.surf(M{1},P);
    [s2,d2] = H.surf(M{2},P);
    fprintf('  %8.1f | %10.3f %12.2f | %10.3f %12.2f%s\n', dhi, s1, d1, s2, d2, ...
            tern(abs(d1-dhi)<1e-6,'   <- A still at bound',''));
end
fprintf('\n  Delta is the delay IN EXCESS of the assumed 5 ms, so a "Delta used" of\n');
fprintf('  1.50 means the model predicts a total neural delay of 6.50 ms.\n');

% ---- 2. shoulder, cell by cell ---------------------------------------------
fprintf('\n===== SHOULDER BY CELL (2nd peak / main; 0 = single-peaked) =====\n');
for a=1:2
    m = M{a};
    fprintf('\n  ARM %s   mean %.4f\n', arms{a,1}, m.shoulder);
    fprintf('  %8s |', 'f (Hz)');
    fprintf(' %7.0f', m.slv); fprintf('   <- level (dB SPL)\n');
    for i=1:numel(m.f)
        fprintf('  %8.0f |', m.f(i));
        for j=1:numel(m.slv), fprintf(' %7.3f', m.sho(i,j)); end
        fprintf('\n');
    end
    ok = isfinite(m.sho);
    fprintf('  cells > 0.5 (detector can flip): %d/%d | > 0.9 (coin flip): %d/%d | max %.3f\n', ...
            sum(m.sho(ok)>0.5), sum(ok(:)), sum(m.sho(ok)>0.9), sum(ok(:)), max(m.sho(ok)));
end

% ---- 3. LIKE-FOR-LIKE: score both arms on the cells BOTH measure cleanly ----
% The cell-by-cell breakdown makes the head-to-head surface comparison suspect.
% Arm B's 0.536 ms is computed from 16 latencies of which 5 sit at shoulder
% 0.74-0.99, i.e. cells where the detector can take either peak; arm A's 1.342 ms
% is computed from 16 clean ones (max 0.023). Comparing them at face value credits
% arm B with an advantage that may be the detector reading the earlier peak, which
% is the exploit the shoulder guard exists to price out.
%
% Restricting BOTH arms to the cells where BOTH are single-peaked removes the
% question. Implemented by blanking the excluded cells to NaN and calling the
% SHIPPED surf_term, which already drops non-finite cells, rather than by
% reimplementing the band loss with a mask argument.
fprintf('\n===== LIKE-FOR-LIKE: only cells where BOTH arms are single-peaked =====\n');
bad = (M{1}.sho > 0.5) | (M{2}.sho > 0.5) | ~isfinite(M{1}.sho) | ~isfinite(M{2}.sho);
fprintf('  excluded %d of %d cells (shoulder > 0.5 in either arm)\n', sum(bad(:)), numel(bad));
P.dlo=0; P.dhi=1.5;
for a=1:2
    mm = M{a}; mm.lat(bad) = NaN;
    [sf,db,nv] = H.surf(mm,P);
    [sa,~,~]   = H.surf(M{a},P);
    fprintf('  ARM %-12s all 16 cells %.3f ms | clean %d cells %.3f ms | Delta %.2f\n', ...
            arms{a,1}, sa, nv, sf, db);
end
fprintf('  (if arm B''s advantage shrinks here, part of it was the detector, not the fit)\n');

save(fullfile(RUNDIR,'dsweep.mat'),'M','bad');
disp('DSWEEP_DONE');
function s=tern(c,a,b), if c, s=a; else, s=b; end, end
