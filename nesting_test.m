% DO THE CHAMBER MODELS NEST?  (SN's principle: a 3-chamber model should never
% be worse than a 2-chamber one, because it can emulate it; if it is worse,
% the APPROACH is at fault.)
%
% STRUCTURAL FINDING FIRST: they do NOT nest as coded. force_cp uses two
% different force laws --
%   m<3 : active force driven by the RELATIVE displacement d3 = d1 - d2
%   m==3: active force driven by d2 ALONE (k_act = gh*k3 - gam*k4)
% and pa.chsz scales only the FLUID coupling, never force_cp. So no chsz
% setting reduces m=3 to m=2. The chamber-count comparison has therefore been
% varying TWO things at once (force law + chamber count).
%
% pa.m3form=1 gives the 3-chamber the m=2-style force law, so the two effects
% can finally be separated:
%   A  m=2 native                                  REFERENCE
%   B  m=3 native            (d2 law, SS open)     current champion structure
%   C  m=3 chsz=[1 0 1]      (d2 law, SS closed)   chamber collapse only
%   D  m=3 chsz=[1 0 1] m3form=1                   TRUE reduction attempt
%   E  m=3 native chsz, m3form=1                   force law only
% If D ~ A, the models nest once the force law is matched and SN's principle
% applies (=> the multi-chamber fits are merely stuck). If D differs from A,
% something else blocks the reduction.
%
% Scored with score26 'fast' -- the common yardstick, so these are directly
% comparable to the champion table.

cfg = { 'A  m=2 native            ', mk(2, [], [])
        'B  m=3 native            ', mk(3, [], [])
        'C  m=3 chsz=[1 0 1]      ', mk(3, [1 0 1], [])
        'D  m=3 chsz=[1 0 1] +m2law', mk(3, [1 0 1], 1)
        'E  m=3 native chsz +m2law', mk(3, [], 1) };

R = struct('tag',{},'S',{});
fprintf('\n  config                     | maperr | range mono | contr | amp    | maxRe\n');
fprintf('%s\n', repmat('-',1,82));
for i = 1:size(cfg,1)
    tag = cfg{i,1}; pa = cfg{i,2};
    try
        S = score26(pa, 'fast', false);
    catch e
        fprintf('  %-26s | FAILED: %s\n', tag, e.message); continue
    end
    fprintf('  %-26s | %6.1f | %5.2f %-4s | %5.1f | %+6.2f | %+7.1f\n', ...
            tag, S.maperr, S.bf_range, S.bf_mono, S.contrast, S.amp_gain, S.maxRe);
    R(end+1).tag = tag; R(end).S = S; %#ok<SAGROW>
end
save('nesting_test.mat','R');
fprintf(['\n  m=2 champion (parfit26_lowdamp) maperr = 191.3 ; m=1 champion = 74.6\n' ...
         '  NESTING HOLDS if D matches A. If D differs, the reduction is blocked by\n' ...
         '  something beyond the force law and the chamber sizes.\n']);
disp('NESTING_TEST_DONE');

function pa = mk(nch, chsz, m3f)
pa = modpar26(nch);
if (~isempty(chsz)), pa.chsz = chsz; end
if (~isempty(m3f)),  pa.m3form = m3f; end
end
