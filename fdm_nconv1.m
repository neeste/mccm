% IS m=1's maperr ADVANTAGE REAL, OR IS IT GRID-SPECIFIC TOO?
%
% WHY THIS MATTERS NOW. SN raises m=1 + an internal d3 as a fallback if the
% multi-chamber models prove intractable. The case for it rests on m=1's maperr
% of 74.6, roughly SEVEN TIMES better than m=3b (499.3) or m=4 (525.2). But that
% 74.6 was measured at n=1401 with the same metric fdm_nconv just showed to be
% U-SHAPED for both multi-chamber models:
%     m=3  956.8 / 401.8 / 438.0 / 499.3 / 636.7 / 790.4 / 1114.9
%     m=4  762.5 / 557.0 / 525.1 / 525.2 / 584.7 / 663.9 /  841.1
% Both have a MINIMUM and worsen in both directions -- no converged value
% exists, and each minimum sits near the grid where that model was fitted.
%
% SO THE QUESTION IS NOT which model scores best at n=1401. It is whether m=1's
% advantage SURVIVES a change of grid. If m=1 is also U-shaped with a minimum at
% 1401, its 74.6 is a fitted-grid artifact of the same kind and the comparison
% is meaningless. If m=1 is FLAT where the multi-chamber models are not, that is
% strong evidence it is both genuinely better AND better conditioned -- which
% would make SN's fallback the serious option rather than the contingency.
%
% m=2 is included because SN notes m=1 already represents 2-chamber physics
% under a symmetry constraint, and m=1 == m=2 was verified exactly earlier in
% this project. They should track each other; if they do not, that identity has
% broken somewhere and would need chasing before anything else here is read.
%
% CHEAP: fdm26 only, no timestep, no corrector. Seconds per point.

NN = [351 701 1051 1401 2101 2801 4201];

fprintf('\n  m=1 reference at n=1401: maperr 74.6 (the best in the project)\n');
fprintf('  m=3 and m=4 are U-SHAPED in n -- minimum near their fitted grid, then\n');
fprintf('  worsening in both directions. Does m=1 do the same?\n\n');
fprintf('  n     | m=1 maperr | m=2 maperr | spread so far\n');
fprintf('%s\n', repmat('-',1,52));
v1 = []; R = struct('n',{},'m1',{},'m2',{});
for nn = NN
    v = nan(1,2);
    for j = 1:2
        pa = modpar26(j);
        if (nn ~= pa.n)
            try, pa = setn(pa, nn); catch, continue; end
        end
        try
            evalc('Rf = fdm26(struct(''pa'',pa));');
            if (isstruct(Rf) && isfield(Rf,'maperr')), v(j) = Rf.maperr; end
        catch
        end
    end
    if (isfinite(v(1))), v1(end+1) = v(1); end %#ok<SAGROW>
    sp = NaN; if (numel(v1)>=2), sp = max(v1)-min(v1); end
    fprintf('  %5d | %10.1f | %10.1f | %12.1f\n', nn, v(1), v(2), sp);
    R(end+1).n=nn; R(end).m1=v(1); R(end).m2=v(2); %#ok<SAGROW>
    save('fdm_nconv1.mat','R');
end
fprintf(['\n  m=1 FLAT while m=3/m=4 are U-shaped -> its advantage is real and it is\n' ...
         '    the better-conditioned model. SN''s m=1 + internal d3 route becomes\n' ...
         '    the serious option, not the fallback.\n' ...
         '  m=1 ALSO U-shaped with a minimum at 1401 -> the 74.6 is a fitted-grid\n' ...
         '    artifact of the same kind, every model in this project shares one\n' ...
         '    disease, and no cross-model comparison made so far can be trusted.\n' ...
         '  m=1 and m=2 DIVERGING from each other -> the verified m=1 == m=2\n' ...
         '    identity has broken, and that must be chased before anything else.\n']);
disp('FDM_NCONV1_DONE');
