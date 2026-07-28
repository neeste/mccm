function R = constr_search(grid)
% CONSTR_SEARCH  Month-2 constrained search over dispersion parameters.
%   Objective : fd group-delay slope d -> 0.413 (ABR target).
%   Constraints: (1) tuning/MAP fit preserved  (maperr <= 1.6x baseline)
%                (2) sub-critical              (coupeig oscillatory max Re < 0)
%   Two stage: cheap fd screen (d + maperr) over the grid, then coupeig on the
%   best tuning-preserving candidates. Reveals the slope-vs-tuning tradeoff and
%   whether any single-param point reaches target while keeping the tuning fit.
%
%   grid.xtap, grid.xtex (etc.) give the sweep values; default xtap x xtex.

if (nargin<1)
    grid.xtap = [0.234 0.28 0.32 0.36 0.40];
    grid.xtex = [4 6 8];
end
base = modpar26(1);
Rb = fdm26(struct('pa',base));
maptol = 1.6*Rb.maperr;
fprintf('baseline: d=%.3f  maperr=%.1f  (tuning tol=%.1f)   target d=0.413\n\n', Rb.d, Rb.maperr, maptol);

% ---- Stage 1: fd screen over the grid ----
xt = grid.xtap; xe = grid.xtex;
rows = zeros(numel(xt)*numel(xe), 5);   % [xtap xtex d maperr tipgain]
i = 0;
for a = xt
    for b = xe
        pa = base; pa.xtap = a; pa.xtex = b;
        Rr = fdm26(struct('pa',pa));
        i = i+1; rows(i,:) = [a b Rr.d Rr.maperr Rr.tipgain];
    end
end
feas = rows(:,4) <= maptol;             % tuning-preserving
obj  = abs(rows(:,3) - 0.413);
[~,ord] = sortrows([~feas obj]);        % feasible first, then closest to target
rows = rows(ord,:); feas = feas(ord); obj = obj(ord);

fprintf('-- fd screen (sorted; * = tuning-preserving, maperr<=%.0f) --\n', maptol);
fprintf('%6s %5s %8s %9s %7s\n','xtap','xtex','slope_d','maperr','tip');
for k = 1:size(rows,1)
    fprintf('%6.3f %5g %8.3f %9.1f %7.1f  %s\n', rows(k,1),rows(k,2),rows(k,3),rows(k,4),rows(k,5), tern(feas(k),'*',' '));
end

% ---- Stage 2: coupeig sub-criticality on the top tuning-preserving candidates ----
cand = find(feas, 3);
fprintf('\n-- coupeig sub-criticality on top %d tuning-preserving candidates --\n', numel(cand));
best = [];
for k = cand'
    pa = base; pa.xtap = rows(k,1); pa.xtex = rows(k,2);
    evalc('S = tdm26(''coupeig'', struct(''pa'',pa));');   % suppress verbose eig dump
    ok = S.maxRe_osc < 0;
    fprintf('xtap=%.3f xtex=%g:  d=%.3f  maperr=%.1f  osc-maxRe=%+.1f  => %s\n', ...
        rows(k,1),rows(k,2),rows(k,3),rows(k,4),S.maxRe_osc, tern(ok,'sub-critical','UNSTABLE'));
    if (ok && isempty(best)), best = [rows(k,1) rows(k,2) rows(k,3) rows(k,4)]; end
end
if (~isempty(best))
    fprintf('\nBEST feasible: xtap=%.3f xtex=%g -> d=%.3f (from %.3f; target 0.413), maperr=%.1f (base %.1f)\n', ...
        best(1),best(2),best(3),Rb.d,best(4),Rb.maperr);
end
R.rows = rows; R.feas = feas; R.baseline = Rb; R.maptol = maptol; R.best = best;
end

function s = tern(c,a,b), if c, s=a; else, s=b; end, end
