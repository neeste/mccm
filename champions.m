% CHAMPION BASELINE: score the current best configuration of each chamber
% count on ONE common yardstick (score26), so future changes -- parameter or
% structural -- can be measured against a fixed reference and either banked or
% rejected on evidence. This is the table the "persistent incremental
% improvement" goal needs before any new fitting starts.
%
% Saves after EACH champion (champions.mat), so an interruption costs at most
% one config -- this session has already lost two runs to restarts.
%
% NOTE ON COMPARABILITY: the historical "bandRMS" figures (0.931 / 0.948 /
% 1.156) came from surf_term, a LOCAL subfunction of parfit26 and therefore not
% callable here. score26 uses the standalone abr_surface_obj instead, so the
% surface numbers below are a NEW consistent baseline -- compare them with each
% other, not with the old bandRMS values.

cfg = {};
% ---- m=1 -------------------------------------------------------------------
p1 = modpar26(1); src1 = 'modpar26(1) native';
if (exist('refit_c1broad.mat','file'))
    L=load('refit_c1broad.mat');
    if (isfield(L,'R') && isfield(L.R,'pa')), p1=L.R.pa; src1='refit_c1broad.mat'; end
end
cfg{end+1} = {1, src1, p1};
% ---- m=2: best surface fit to date ----------------------------------------
p2 = modpar26(2); src2 = 'modpar26(2) native';
if (exist('parfit26_lowdamp.mat','file'))
    L=load('parfit26_lowdamp.mat');
    if (isfield(L,'R') && isfield(L.R,'pa')), p2=L.R.pa; src2='parfit26_lowdamp.mat'; end
end
cfg{end+1} = {2, src2, p2};
% ---- m=3: the (corrected) milestone fit ------------------------------------
p3 = modpar26(3); src3 = 'modpar26(3) native';
if (exist('refit_c3_map.mat','file'))
    L=load('refit_c3_map.mat');
    if (isfield(L,'R') && isfield(L.R,'pa')), p3=L.R.pa; src3='refit_c3_map.mat'; end
end
cfg{end+1} = {3, src3, p3};
% ---- m=4: best REAL config = unfit m2o x32 --------------------------------
% (parfit26_c4map did not decompress the map; parfit26_c4map_m2 reached a good
%  fdm26 maperr but DIVERGES in the time domain, so neither is usable.)
p4 = modpar26(4); p4.m2o = modpar26(4).m2o * 32;
cfg{end+1} = {4, 'modpar26(4) with m2o x32', p4};

C = struct('nch',{},'src',{},'S',{});
for i = 1:numel(cfg)
    nch = cfg{i}{1}; src = cfg{i}{2}; pa = cfg{i}{3};
    fprintf('\n############ CHAMPION m=%d  <- %s ############\n', nch, src);
    t0 = tic;
    try
        S = score26(pa, 'full');
    catch e
        fprintf('  score26 FAILED: %s\n', e.message);
        S = struct('m',nch,'err',e.message);
    end
    S.src = src; S.wall = toc(t0);
    fprintf('  (%.0f s)\n', S.wall);
    C(end+1).nch = nch; C(end).src = src; C(end).S = S; %#ok<SAGROW>
    save('champions.mat','C');            % bank after every config
end

% ---- summary table ---------------------------------------------------------
fprintf('\n\n================= CHAMPION SUMMARY =================\n');
fprintf(' m | maperr | maxRe(all) | range(oct) mono | contr | amp   | surf-shape |  d    |  b    | level%%/dB | OAE lat  ratio\n');
for i = 1:numel(C)
    S = C(i).S;
    f = @(fn,fmt) fmtf(S,fn,fmt);
    fprintf(' %d | %s | %s | %s %-4s | %s | %s | %s | %s | %s | %s | %s %s\n', ...
        C(i).nch, f('maperr','%6.1f'), f('maxRe','%+10.1f'), f('bf_range','%5.2f'), ...
        gets(S,'bf_mono'), f('contrast','%5.1f'), f('amp_gain','%+5.2f'), ...
        f('surf_resid','%10.3f'), f('d','%5.3f'), f('b','%5.2f'), ...
        f('level_pct','%9.2f'), f('oae_lat','%7.2f'), f('oae_ratio','%5.2f'));
end
fprintf('\n targets:            d 0.39-0.41 | b 11.99-13.27 | level 1.62 %%/dB | OAE ratio ~2\n');
fprintf(' m=2 native CF range is ~5.9 oct -- the tonotopic span a real cochlea needs.\n');
disp('CHAMPIONS_DONE');

function s = fmtf(S,fn,fmt)
if (isfield(S,fn) && ~isempty(S.(fn)) && isnumeric(S.(fn)) && isfinite(S.(fn)(1)))
    s = sprintf(fmt, S.(fn)(1));
else
    s = sprintf(sprintf('%%%ds', numel(sprintf(fmt,0))), '--');
end
end
function v = gets(S,fn)
if (isfield(S,fn) && ischar(S.(fn))), v = S.(fn); else, v = '?'; end
end
