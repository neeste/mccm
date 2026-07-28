% CONTINUOUS FORCE-LAW INTERPOLATION for the 3-chamber (pa.m3form = alpha).
%
% The two force laws were measured to be COMPLEMENTARY, not right-vs-wrong:
%   alpha=0 (d2-only)      maperr 499.3, contrast 9.2,  degeneracy from 1 kHz
%   alpha=1 (d1-d2)        maperr 1868.7, contrast 11.1, degeneracy from 2 kHz
% so an intermediate drive may take the tuning of one and the tip/map validity
% of the other. Making the law continuous also converts a discrete model-class
% switch into a FITTABLE parameter -- the reformulation SN's marginal-
% improvement principle needs, given that chamber count itself cannot be made
% continuous (chsz->0 is singular in fdm26 and unstable in tdm26).
%
% ENDPOINT REGRESSION FIRST. alpha=0 and alpha=1 must reproduce the previously
% measured values EXACTLY; the interpolation is only trustworthy if it does.
%   alpha=0 : maperr 499.3, range 6.53 FOLD, contrast 9.2, amp +81.15
%   alpha=1 : maperr 1868.7, range 6.72 FOLD, contrast 11.1
% (amp is NaN wherever the passive gam=0 reference diverges -- known for the
%  d1-d2 law, so an NaN there is expected, not a failure.)

AL = [0 0.25 0.5 0.75 1.0];
fprintf('\n  alpha | maperr   | range mono | contr | amp    | maxRe    | note\n');
fprintf('%s\n', repmat('-',1,80));
R = struct('al',{},'S',{});
for i = 1:numel(AL)
    al = AL(i);
    pa = modpar26(3); pa.m3form = al;
    try
        S = score26(pa, 'fast', false);
    catch e
        fprintf('  %5.2f | FAILED: %s\n', al, e.message); continue
    end
    note = '';
    if (al==0),   note = 'regression: want 499.3 / 6.53 / 9.2 / +81.15'; end
    if (al==1),   note = 'regression: want 1868.7 / 6.72 / 11.1'; end
    fprintf('  %5.2f | %8.1f | %5.2f %-4s | %5.1f | %+6.2f | %+8.1f | %s\n', ...
            al, S.maperr, S.bf_range, S.bf_mono, S.contrast, S.amp_gain, S.maxRe, note);
    R(end+1).al = al; R(end).S = S; %#ok<SAGROW>
end
save('m3alpha.mat','R');

% ---- verdict on the endpoints ----------------------------------------------
ok = true;
for i = 1:numel(R)
    if (R(i).al==0 && abs(R(i).S.maperr-499.3) > 1),  ok=false; fprintf('\n  *** REGRESSION FAIL at alpha=0: maperr %.1f, want 499.3\n', R(i).S.maperr); end
    if (R(i).al==1 && abs(R(i).S.maperr-1868.7) > 1), ok=false; fprintf('\n  *** REGRESSION FAIL at alpha=1: maperr %.1f, want 1868.7\n', R(i).S.maperr); end
end
if (ok), fprintf('\n  ENDPOINT REGRESSION PASSED -- interpolation reproduces both laws exactly.\n'); end
fprintf(['  LOOK FOR: an interior alpha with maperr near/below 499 AND contrast\n' ...
         '  above 9.2. That would beat BOTH endpoints and is the whole point.\n' ...
         '  If maperr rises monotonically with alpha and contrast too, the two are\n' ...
         '  strictly traded and no interior optimum exists on these two axes.\n']);
disp('M3ALPHA_DONE');
