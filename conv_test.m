% GRID CONVERGENCE: is n=701 as good as n=1401, and how much time does it save?
%
% SN selected n=1401 to maintain numerical stability when the tuning becomes
% very sharp. Two facts established before running this:
%  (a) pa.dt = 2e-6 is FIXED in modpar26.m and does NOT scale with n, so halving
%      n doubles dx while dt stays put. That moves the explicit scheme FURTHER
%      from its CFL limit, not closer. The risk is therefore SPATIAL RESOLUTION
%      of a narrow near-CF envelope, not time-step stability.
%  (b) score26 previously hardcoded n=1401 place indices; they are now computed
%      as fractions of pa.n, so the SAME FRACTIONAL PLACES are compared at both
%      grid sizes. Without that fix any n=701 number would be meaningless.
%  (c) pa carries n-LENGTH PROFILE ARRAYS (gampro, synpro) in both saved fits
%      and a fresh modpar26(k). Setting pa.n alone leaves them at the old length
%      and tdm26 throws immediately. setn.m resamples them (and rescales isv).
%      The first run of this test failed entirely for that reason.
%
% Sharply tuned configs are the risk cases, so m=3 native (+81 dB, the sharpest
% amplifier measured) is the strongest test of SN's stated concern. m=4 at
% +2.46 dB is NOT sharply tuned, so agreement there proves little on its own.

L = load('refit_c1broad.mat');
b2 = modpar26(4).m2o;
p4 = modpar26(4); p4.m2o = b2*32;
cfg = { 'm=1 champion  (amp +50.7)', L.R.pa
        'm=3 native    (amp +81.2)', modpar26(3)
        'm=4 m2o x32   (amp  +2.5)', p4 };

NN = [1401 701];
fprintf('\n  config                    |    n | maperr | range mono fold | contr | amp    | maxRe    | sec\n');
fprintf('%s\n', repmat('-',1,104));
R = struct('tag',{},'n',{},'S',{},'sec',{});
for i = 1:size(cfg,1)
    for j = 1:numel(NN)
        pa = setn(cfg{i,2}, NN(j));   % resamples gampro/synpro/isv; setting
                                      % pa.n alone throws a size mismatch
        t0 = tic;
        try
            S = score26(pa, 'fast', false);
        catch e
            fprintf('  %-25s | %4d | FAILED: %s\n', cfg{i,1}, NN(j), e.message); continue
        end
        sec = toc(t0);
        fprintf('  %-25s | %4d | %6.1f | %5.2f %-4s %.2f | %5.1f | %+6.2f | %+8.1f | %4.0f\n', ...
                cfg{i,1}, NN(j), S.maperr, S.bf_range, S.bf_mono, S.bf_fold, ...
                S.contrast, S.amp_gain, S.maxRe, sec);
        R(end+1).tag=cfg{i,1}; R(end).n=NN(j); R(end).S=S; R(end).sec=sec; %#ok<SAGROW>
        save('conv_test.mat','R');
    end
end

% ---- differences and speedup ----------------------------------------------
fprintf('\n  config                    | d(maperr) | d(range) | d(contr) | d(amp) | speedup\n');
fprintf('%s\n', repmat('-',1,84));
for i = 1:size(cfg,1)
    k = find(strcmp({R.tag}, cfg{i,1}));
    if (numel(k) < 2), continue; end
    a = R(k(1)).S; b = R(k(2)).S;          % a = 1401, b = 701
    fprintf('  %-25s | %+9.1f | %+8.2f | %+8.1f | %+6.2f | %5.2fx\n', ...
            cfg{i,1}, b.maperr-a.maperr, b.bf_range-a.bf_range, ...
            b.contrast-a.contrast, b.amp_gain-a.amp_gain, R(k(1)).sec/R(k(2)).sec);
end
fprintf(['\n  ADOPT n=701 for exploration only if the differences are small on the\n' ...
         '  SHARPLY TUNED config (m=3, +81 dB). Agreement on m=4 (+2.5 dB) does not\n' ...
         '  license n=701 for the high-gain regime we are trying to reach, since\n' ...
         '  that is exactly where the near-CF envelope narrows.\n']);
disp('CONV_TEST_DONE');
