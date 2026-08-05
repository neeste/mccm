% M4DIAG -- WHY does an m=4 fit fail the getpar/setpar round-trip by 136.5?
%
% Established 2026-08-04: the discrepancy is PRE-EXISTING, not caused by the
% k5/r5/m5 append (it survives with the m>=4 gate applied, so that path is
% byte-identical to before). What has NOT been established is the cause, and the
% standing hypothesis -- "k5q/r5q/m5q are not transported" -- is a guess.
%
% The mechanism is structural: setpar_l REBUILDS pa from base=modpar26(m) and
% writes back only what pv carries (30 impedance + chsz + gpgrade + hbsc, plus
% the six third-DOF entries at m<4). ANY field of a saved fit that differs from
% its modpar26 default and is not in that list is SILENTLY RESET on every
% evaluation. So the question is answerable by comparison, with no model runs:
% list the fields that differ, subtract the ones pv transports, and what remains
% is the complete candidate set. Then confirm by restoring them one at a time --
% a cause that does not reproduce the 136.5 is not the cause.
PROJ = '/Users/neely/Library/CloudStorage/OneDrive-FatherFlanagan''sBoysHome/STN/mod/mccm';
addpath(PROJ); cd(PROJ);
FIT = 'refit_m4_full.mat';

L=load(FIT); pa=L.R.pa; if (~isfield(pa,'m')), pa=L.R; end
base = modpar26(pa.m);
H = parfit26('handles');
nc = numel(pa.chsz);
carried = [H.parnames(), {'chsz','gpgrade','gampro','hbsc','m'}];

fprintf('\n===== M4DIAG: %s (m=%d) =====\n', FIT, pa.m);

% ---- 1. every field that differs from the modpar26 default ----
fn = fieldnames(pa); dif = {};
for i=1:numel(fn)
    f = fn{i};
    if (~isfield(base,f)), dif{end+1}=f; continue; end %#ok<SAGROW>
    a=pa.(f); b=base.(f);
    if (~isnumeric(a) && ~islogical(a)), continue; end
    if (~isequal(size(a),size(b))), dif{end+1}=f; continue; end %#ok<SAGROW>
    if (any(a(:)~=b(:))), dif{end+1}=f; end %#ok<SAGROW>
end
lost = setdiff(dif, carried);
fprintf('\n  fields differing from modpar26(%d): %d\n', pa.m, numel(dif));
fprintf('  of those, NOT transported by pv (silently reset every eval): %d\n', numel(lost));
for i=1:numel(lost)
    f=lost{i}; a=pa.(f);
    if (isfield(base,f)), b=base.(f); else, b=NaN; end
    if (isscalar(a) && isscalar(b))
        fprintf('    %-10s saved %-14.6g default %-14.6g\n', f, a, b);
    else
        fprintf('    %-10s saved [%s] default [%s]\n', f, num2str(size(a)), ...
                tern(isfield(base,f),num2str(size(b)),'ABSENT'));
    end
end
if (isempty(lost))
    fprintf('  NONE -- the round-trip loss is not a dropped field. Look elsewhere.\n');
end

% ---- 2. the damage, and whether restoring the candidates repairs it ----
pe = H.setpar(base, H.getpar(pa), []);
e0 = fdm26(struct('pa',pa)); e0=e0.maperr;
e1 = fdm26(struct('pa',pe)); e1=e1.maperr;
fprintf('\n  saved maperr %.4f | after round-trip %.4f | delta %.4f\n', e0, e1, abs(e0-e1));

fprintf('\n  restoring each lost field ALONE onto the round-tripped model:\n');
fprintf('  %-12s %12s %12s\n','field','maperr','remaining');
for i=1:numel(lost)
    q = pe; q.(lost{i}) = pa.(lost{i});
    try
        r = fdm26(struct('pa',q));
        fprintf('  %-12s %12.4f %12.4f%s\n', lost{i}, r.maperr, abs(e0-r.maperr), ...
                tern(abs(e0-r.maperr)<1e-6,'   <== THIS IS THE CAUSE',''));
    catch e, fprintf('  %-12s FAILED %s\n', lost{i}, e.message(1:min(40,end))); end
end

% ---- 3. all of them together: is the set COMPLETE? ----
q = pe; for i=1:numel(lost), q.(lost{i}) = pa.(lost{i}); end
r = fdm26(struct('pa',q));
fprintf('\n  ALL restored: maperr %.4f | remaining delta %.4e %s\n', r.maperr, abs(e0-r.maperr), ...
        tern(abs(e0-r.maperr)<1e-9,'  <== set is COMPLETE', ...
                                    '  <== set is INCOMPLETE, something else is lost too'));
save('m4diag.mat','lost','dif','e0','e1');
disp('M4DIAG_DONE');
function s=tern(c,a,b), if c, s=a; else, s=b; end, end
