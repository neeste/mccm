% IS THE FORWARD-LATENCY SLOPE CONVERGED IN n?
%
% WHY. eval_cost found m=4's slope reads 0.333 at n=701 and 0.157 at n=1401 --
% it MORE THAN DOUBLES between grids. A quantity that moves that far is not
% converged, and there is no reason to assume 1401 is the resting place. This
% undercuts the conclusion that m=4's slope is "badly broken" at 0.157 against
% the 0.413 target: the n=701 value is much CLOSER to target than the n=1401
% one, so the direction of convergence decides whether m=4 has a real latency
% defect or a discretization artifact.
%
% THE CONTROL IS THE POINT. m=3b is swept at all three grids too, and it is
% cheap (abr ~100 s at n=1401 against m=4's 2290 s). If m=3b's 0.400 is stable
% across grids, the problem is specific to the nested m=4 build. If m=3b ALSO
% moves, then slope is grid-sensitive generally, every slope in this project is
% suspect -- including the 0.400 that appears to sit on target and the 0.413
% target itself, which was established from model output -- and that is a much
% larger problem than m=4.
%
% KNOWN (eval_cost): m=4 n=701 slope 0.333 maperr 557.0
%                    m=4 n=1401 slope 0.157 maperr 525.2
% Those are re-measured here rather than assumed, so the whole table comes from
% one run under one set of conditions.
%
% STABILITY RISK AT n=2801. The time march uses a FIXED dt=2e-6. Halving dx
% without touching dt has previously broken this integrator -- stiffening k1o
% sent every run non-finite for the same reason. So a divergent n=2801 row is a
% possible and INFORMATIVE outcome: it would mean the convergence question
% cannot be answered by refining alone, and the fixed dt would have to be
% addressed first. Guarded and reported rather than left to throw.
%
% setn() is used for every n change -- it resamples gampro/synpro and rescales
% isv. Setting pa.n directly throws or silently corrupts.
%
% COST m=3b is cheap at all three grids; m=4 at n=2801 is the expensive row,
% plausibly 100-150 min alone. Other jobs are running, so wall time is inflated.

NN = [701 1401 2801];

fprintf('\n  target slope 0.413 | m=3b n=1401 reads 0.400 | m=4 n=1401 reads 0.157\n');
fprintf('  CONTROL FIRST: if m=3b moves across grids, every slope in the project\n');
fprintf('  is suspect, not just m=4''s.\n\n');
fprintf('  model | n    | slope | level_c | maperr | finite\n');
fprintf('%s\n', repmat('-',1,56));
R = struct('m',{},'n',{},'slope',{},'lc',{},'maperr',{});
for nch = [3 4]
    for nn = NN
        pa = modpar26(nch);
        if (nn ~= pa.n)
            try, pa = setn(pa, nn); catch e
                fprintf('  m=%d   | %4d | setn FAILED: %s\n', nch, nn, e.message); continue
            end
        end
        sl = NaN; lc = NaN; mp = NaN; fin = 'n/a';
        try
            evalc('m = abr_metric(pa,false);'); sl = m.slope; lc = m.level_c;
            fin = 'yes'; if (~isfinite(sl)), fin = 'NO'; end
        catch e
            fprintf('  m=%d   | %4d | abr FAILED: %s\n', nch, nn, e.message); fin='THREW';
        end
        try
            S = score26(pa,'fast',false); mp = S.maperr;
        catch
        end
        fprintf('  m=%d   | %4d | %5.3f | %7.2f | %6.1f | %s\n', nch, nn, sl, lc, mp, fin);
        R(end+1).m=nch; R(end).n=nn; R(end).slope=sl; R(end).lc=lc; R(end).maperr=mp; %#ok<SAGROW>
        save('conv_slope.mat','R');
    end
end
fprintf(['\n  READ THE m=3 ROWS FIRST.\n' ...
         '  m=3 stable, m=4 moving -> the nested build has a genuine grid\n' ...
         '    sensitivity, and m=4''s slope must be quoted at a converged n before\n' ...
         '    being called broken.\n' ...
         '  BOTH moving -> slope is grid-sensitive generally and every slope\n' ...
         '    figure in this project needs re-establishing, target included.\n' ...
         '  m=4 settling toward 0.157 as n rises -> the defect is real after all.\n' ...
         '  m=4 settling toward ~0.33 -> much closer to m=3b and to target, and\n' ...
         '    the "badly broken" reading was a discretization artifact.\n']);
disp('CONV_SLOPE_DONE');
