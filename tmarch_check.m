% TIME-MARCH STABILITY CHECK ON THE ohcgain=0.01 CANDIDATE.
%
% The backoff sweep found one configuration that combines a clean 4-chamber map
% with a real amplifier: ohcsgn=-1, ohcgain=0.01, m2o x32, r1 x0.5 gives
% amp +29.49 dB, range 4.32 oct, fold 0.00, contrast 54.4, maxRe +185.6.
% My script called it UNSTABLE, but that verdict rests on a threshold of 24 that
% I derived from only FOUR champions (+0.0 to +23.6). +185.6 sits ~300x BELOW
% the genuinely divergent cases (+62747 and up), and the click ran to
% completion, so the eigenvalue threshold alone cannot settle it.
%
% A short click can miss slow growth. This drives the model with TONE BURSTS at
% several frequencies and levels, which run far longer, and measures whether the
% response GROWS. A self-oscillating model shows late-time growth even when the
% early response looks normal.
%
% TESTS
%  1 finiteness and peak magnitude across a 4x4 frequency/level grid (tbabr)
%  2 late-vs-early energy ratio of the WNR at each condition: >1 and rising with
%    time indicates growth rather than a decaying response
%  3 the ohcgain=0.005 neighbour (maxRe +0.0, amp -0.30) as a STABLE CONTROL,
%    so the growth measure is calibrated against a model known to be quiet

b2 = modpar26(4).m2o; br1 = modpar26(4).r1o;
mk = @(og) setfield(setfield(setfield(setfield(modpar26(4), ...
        'm2o', b2*32), 'r1o', br1*0.5), 'ohcsgn', -1), 'ohcgain', og); %#ok<SFLD>

for og = [0.01 0.005]
    pa = mk(og); pa.hbmode = 'bm';
    fprintf('\n===== ohcgain = %g =====\n', og);
    % ---- 1: full tone-burst grid, finiteness and latency sanity -------------
    try
        evalc('T = tdm26(''tbabr'', pa, 0, 0);');
        nfin = nnz(~isfinite(T.lat));
        fprintf('  tbabr grid: %d/%d cells non-finite latency\n', nfin, numel(T.lat));
        fprintf('  latency (ms) by freq x level:\n');
        for j = 1:size(T.lat,1)
            fprintf('    %5.2f kHz |', T.f(j));
            fprintf(' %7.2f', T.lat(j,:));
            fprintf('\n');
        end
    catch e
        fprintf('  tbabr FAILED: %s\n', e.message);
    end
    % ---- 2: late-vs-early growth in the WNR at several conditions -----------
    fprintf('  growth check (late/early rms of WNR; ~1 = steady, >>1 = growing)\n');
    for fr = [1 2 4]
        for lv = [40 80]
            p.fr = fr; p.lv = lv; p.pa = pa;
            try
                evalc('S = tdm26(''wnr1'', p, 0, 0);');
                w = S.wnr(:); w = w(isfinite(w));
                nq = max(4, floor(numel(w)/4));
                e1 = sqrt(mean(w(1:nq).^2));            % first quarter
                e4 = sqrt(mean(w(end-nq+1:end).^2));    % last quarter
                fprintf('    %g kHz %2d dB : late/early = %8.3g   max|WNR| %9.3e\n', ...
                        fr, lv, e4/max(e1,eps), max(abs(w)));
            catch e
                fprintf('    %g kHz %2d dB : FAILED %s\n', fr, lv, e.message);
            end
        end
    end
end
fprintf(['\n  READ: ohcgain=0.005 is the quiet control (maxRe +0.0). If 0.01 shows\n' ...
         '  a similar late/early ratio and finite latencies everywhere, maxRe +185.6\n' ...
         '  is not disqualifying and my threshold of 24 was too strict. If 0.01\n' ...
         '  grows late while 0.005 does not, it is self-oscillating and unusable.\n']);
disp('TMARCH_CHECK_DONE');
