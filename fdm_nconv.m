% IS maperr CONVERGED IN n? (frequency domain only)
%
% WHY THIS IS SEPARATE FROM EVERYTHING ELSE TODAY. maperr comes from fdm26, a
% FREQUENCY-DOMAIN solver: no timestep, no corrector, no nimp. nimp_conv proved
% it -- maperr read exactly 499.3 for nimp 1, 2, 3 and 5. So its variation with n
% is a PURE SPATIAL effect, unrelated to dt/tau or to the corrector, and the
% time-domain diagnosis does not touch it.
%
% WHAT conv_slope FOUND. maperr ran 401.8 -> 499.3 -> 790.4 for m=3 across
% n = 701/1401/2801: rising by 24% then 58%, monotonically, AWAY from the
% targets. A discretization error that grows with resolution is not a
% discretization error.
%
% MY FIRST EXPLANATION WAS WRONG. I guessed maperr was an accumulated sum over
% places and so scaled with n. fdm26.m:598 shows otherwise: it is evaluated over
% SEVEN FIXED frequencies, flmap = 500*2.^(-1:5) = 250 Hz to 16 kHz. The count
% does not depend on n, so this is not an accumulation artifact -- the model's
% tuning at those seven frequencies genuinely moves with spatial resolution.
%
% THE ODD PART, and the reason this matters. n=701 scores BETTER (401.8) than
% the n=1401 (499.3) at which the model was fitted. A metric the parameters were
% tuned against should not improve when the grid is coarsened. Either the
% solution is far from converged, or maperr is measuring something that does not
% behave monotonically in resolution.
%
% THIS IS CHEAP. fdm26 alone is seconds to a couple of minutes, against 38
% MINUTES for a single m=4 abr_metric. A wide sweep costs less than one m=4
% latency point, so there is no reason not to bracket the question properly.
%
% WHAT TO LOOK FOR
%   maperr FLATTENS at high n  -> there is a converged answer; n=1401 is simply
%                                 too coarse, and the fitted parameters are
%                                 grid-specific
%   maperr KEEPS RISING        -> the frequency-domain solution does not
%                                 converge, and every maperr in this project is
%                                 a property of its grid rather than its model

NN = [351 701 1051 1401 2101 2801 4201];

fprintf('\n  maperr over 7 FIXED frequencies (250 Hz - 16 kHz), so the count does\n');
fprintf('  not scale with n. Known: m=3 gave 401.8 / 499.3 / 790.4 at 701/1401/2801.\n');
fprintf('  fdm26 only -- no timestep, no corrector, no nimp involved.\n\n');
fprintf('  n     | m=3 maperr | m=4 maperr | m=3 s | m=4 s\n');
fprintf('%s\n', repmat('-',1,56));
R = struct('n',{},'m3',{},'m4',{});
for nn = NN
    v = nan(1,2); tt = nan(1,2);
    for j = 1:2
        nch = j + 2;                      % 3 then 4
        pa = modpar26(nch);
        if (nn ~= pa.n)
            try, pa = setn(pa, nn); catch, continue; end
        end
        t0 = tic;
        try
            evalc('Rf = fdm26(struct(''pa'',pa));');
            if (isstruct(Rf) && isfield(Rf,'maperr')), v(j) = Rf.maperr; end
        catch
        end
        tt(j) = toc(t0);
    end
    fprintf('  %5d | %10.1f | %10.1f | %5.1f | %5.1f\n', nn, v(1), v(2), tt(1), tt(2));
    R(end+1).n = nn; R(end).m3 = v(1); R(end).m4 = v(2); %#ok<SAGROW>
    save('fdm_nconv.mat','R');
end
fprintf(['\n  If maperr flattens, read off the converged value and note that every\n' ...
         '  comparison made today was at n=1401, i.e. off the converged answer.\n' ...
         '  If it keeps rising, maperr is a property of the grid and not of the\n' ...
         '  model, and the whole map-fitting objective needs reconsidering before\n' ...
         '  any fit is run against it.\n' ...
         '  Either way the RANKINGS taken at fixed n=1401 stay internally valid;\n' ...
         '  it is the absolute values and the cross-n comparisons that are at risk.\n']);
disp('FDM_NCONV_DONE');
