% RLPROBE -- two questions about the same object, the pressure-selector Df.
%
% Q1 (SN): is there value in an RL-WITH-FLUID-FORCE case at m<3?
%     A fluid force on a membrane is a row of Df, i.e. a PRESSURE DIFFERENCE
%     between two chambers. So the question is not whether it is implemented,
%     it is whether it is EXPRESSIBLE at m<3. Part 1 prints the actual Df for
%     each configuration rather than arguing from the source.
%
% Q2: is rlsplit fittable? It is the last unexplored structural candidate
%     (maperr 517.10 unfitted vs 144.02) but it carries a RECORDED HISTORY OF
%     BEING ILL-POSED: divergence at sample 26, invariant under a 40x change in
%     CL area, which is the signature of over-determination rather than of a
%     bad parameter. macro26 has since derived a2 from Df/Dq instead of holding
%     a second copy of the topology, which is the stated fix. Whether it TOOK is
%     an empirical question and it is the go/no-go for spending hours here.
PROJ = '/Users/neely/Library/CloudStorage/OneDrive-FatherFlanagan''sBoysHome/STN/mod/mccm';
addpath(PROJ); cd(PROJ);

fprintf('\n===== PART 1: is a fluid force on the RL EXPRESSIBLE at m<3? =====\n');
CFG = { 'm=1 dof=3',        1, {'dof',3}; ...
        'm=2 dof=3',        2, {'dof',3}; ...
        'm=3 d3int=1',      3, {'d3int',1}; ...
        'm=3 rlsplit',      3, {'d3int',1,'rlsplit',1}; ...
        'm=4',              4, {} };
for i=1:size(CFG,1)
    pa = modpar26(CFG{i,2}); kv = CFG{i,3};
    for j=1:2:numel(kv), pa.(kv{j}) = kv{j+1}; end
    try
        cp = imped26(pa); C = macro_couple(pa, cp);
        fprintf('\n  %-14s nch=%d ndof=%d  dref=[%s]\n', CFG{i,1}, C.nch, C.ndof, ...
                strtrim(sprintf('%d ', C.dref)));
        for k=1:size(C.Df,1)
            nm = {'d1 BM','d2 TM/shear','d3 RL/OC','d4 vent'};
            fprintf('     Df(%d,:) %-11s [%s]%s\n', k, nm{min(k,4)}, ...
                    strtrim(sprintf('%6.2f', C.Df(k,:))), ...
                    tern(k>=2 && all(C.Df(k,:)==0), '   <-- NO fluid force, internal DOF', ''));
        end
    catch e, fprintf('  %-14s FAILED: %s\n', CFG{i,1}, e.message(1:min(60,end))); end
end
fprintf('\n  A membrane needs DISTINCT fluid compartments on its two sides for a\n');
fprintf('  pressure difference to exist. Count the chambers against the rows above.\n');

fprintf('\n===== PART 2: is rlsplit WELL-POSED now? (the go/no-go) =====\n');
L=load('fit_nch3_surface.mat'); pa0=L.R.pa;
if (isfield(pa0,'hbmode')), pa0=rmfield(pa0,'hbmode'); end
V = { 'baseline d3int=1',       struct('d3int',1); ...
      'rlsplit WITHOUT d3int',  struct('rlsplit',1); ...
      'rlsplit + d3int',        struct('d3int',1,'rlsplit',1) };
for i=1:size(V,1)
    pa = pa0; f=fieldnames(V{i,2});
    pa.d3int = 0; pa.rlsplit = 0;
    for j=1:numel(f), pa.(f{j}) = V{i,2}.(f{j}); end
    fprintf('\n  --- %s ---\n', V{i,1});
    try
        cp=imped26(pa); C=macro_couple(pa,cp);
        fprintf('    ndof %d | Df rows %d %s\n', C.ndof, size(C.Df,1), ...
                tern(size(C.Df,1)~=C.ndof,'   <-- MISMATCH: Df has a row the solver will IGNORE',''));
        R=fdm26(struct('pa',pa));
        fprintf('    maperr %.2f\n', R.maperr);
    catch e, fprintf('    fdm26 FAILED: %s\n', e.message(1:min(70,end))); continue; end
    try   % THE ill-posedness test: the old failure was divergence at sample 26
        S=tdm26('wnr1',struct('fr',2,'lv',60,'pa',pa),0,0);
        fprintf('    tdm wnr1 OK: tpk %.2f ms | d1mx %.3e | finite %d\n', ...
                S.tpk, S.dgn.d1mx, isfinite(S.dgn.d1mx));
    catch e, fprintf('    tdm26 FAILED: %s   <-- STILL ILL-POSED\n', e.message(1:min(60,end))); end
    try
        evalc('C2=tdm26(''coupeig'',struct(''pa'',pa));'); s=score26(pa,'fast',false);
        fprintf('    maxRe %+.2f | osc %+.2f | amp %+.2f dB\n', C2.maxRe, C2.maxRe_osc, s.amp_gain);
    catch e, fprintf('    stability FAILED: %s\n', e.message(1:min(60,end))); end
end

fprintf('\n===== PART 3: is rlsplit FITTABLE? =====\n');
H=parfit26('handles');
pa=pa0; pa.d3int=1; pa.rlsplit=1; nc=numel(pa.chsz);
base=modpar26(3); base.d3int=1; base.rlsplit=1;
pe=H.setpar(base,H.getpar(pa),[]);
lost=H.lost(pa,pe,[]);
fprintf('  warm start survives setpar_l: %s\n', tern(isempty(lost),'YES', ...
        ['NO -- lost ' strjoin(lost,', ') ' (pass via opts.pin)']));
try
    e0=fdm26(struct('pa',pe)); e0=e0.maperr;
    pvk=H.getpar(pa); pvk(30+nc+3)=pvk(30+nc+3)*2;
    q=H.setpar(base,pvk,[]); rq=fdm26(struct('pa',q));
    fprintf('  k5o x2 through pv: maperr %.2f -> %.2f (delta %+.3e)%s\n', ...
            e0, rq.maperr, rq.maperr-e0, tern(abs(rq.maperr-e0)<1e-9,'  <-- INERT',''));
catch e, fprintf('  liveness FAILED: %s\n', e.message(1:min(60,end))); end
disp('RLPROBE_DONE');
function s=tern(c,a,b), if c, s=a; else, s=b; end, end
