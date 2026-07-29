function y = hbmix(hb, x1, x2, x3)
%HBMIX  Hair-bundle coordinate from the hb mixing matrix. ONE definition.
%
%   y = HBMIX(hb, d1, d2)      2-column form (no third DOF)
%   y = HBMIX(hb, d1, d2, d3)  3-column form (bundle may reference DOF 3)
%
% Pass hb ALREADY row-selected if you only want some places, e.g.
% hbmix(cp.hb(sv,:), d1, d2, dc).
%
% WHY THIS EXISTS. cp.hb had FIVE live consumers (micro26's MET gain, tdm26's
% IHC hair-bundle velocity, ecochg3, gamcmp, plus construction in macro26).
% Extending it to three columns would have left every consumer silently
% computing the bundle from two of three terms -- wrong only when the third
% column is nonzero, i.e. exactly in the new configuration under test. That is
% the same failure mode as the has3 predicate duplicated at four sites, which
% produced "Unrecognized function or variable 'dc'" when only one was updated.
% Keep this the single definition.
%
% WHAT THE COLUMNS MEAN (set in macro26.m):
%   m<3    [gh  -1   0]   bundle = gh*d1 - d2      d2 is the TM
%   m>=3   [ 0   1   0]   bundle = d2              d2 IS the bundle, TM gone
%   hbrl   [ 0  -1  gh]   bundle = gh*d3 - d2      d3 is the RL (Liu & Neely
%                         2010: MET is driven by RL motion, no explicit TM)
%
% The third column is ZERO in both legacy forms, so this is behaviour-preserving
% until pa.hbrl is set.
if (nargin < 4 || isempty(x3) || size(hb,2) < 3)
    y = hb(:,1) .* x1 + hb(:,2) .* x2;
else
    y = hb(:,1) .* x1 + hb(:,2) .* x2 + hb(:,3) .* x3;
end
end
