% NULLDIR -- move the latency slope WITHOUT paying in level dependence, by using
% a COMBINATION of parameters instead of one at a time.
%
% WHY THE SINGLE-PARAMETER ATLAS ASKED THE WRONG QUESTION. latatlas ranked every
% parameter by selectivity S = |dslope|/|dlevel_c| and found NONE above 0.48:
% every lever buys slope and pays in level dependence. That looks like a wall.
% It is not. Two levers with DIFFERENT trade ratios can be combined so their
% level_c effects cancel while their slope effects do not:
%
%     k3e   dslope +0.25047   dlvlc  -4.0827    ratio -0.0613
%     r2e   dslope +0.15676   dlvlc -10.4040    ratio -0.0151
%
% Set b = -0.3924a and the level_c changes cancel exactly, leaving
% dslope = 0.189a. So a selective direction EXISTS even though no selective
% AXIS does. The single-parameter frame could never have seen it.
%
% This builds the best such direction from the whole atlas: project the slope
% gradient onto the null space of the level_c gradient,
%
%     u = g_slope - (g_slope . g_lvlc / |g_lvlc|^2) g_lvlc
%
% in RELATIVE-perturbation coordinates (the atlas perturbed every parameter by
% the same 5% relative step, so its rows are already commensurately scaled).
%
% LINEARIZATION IS A HYPOTHESIS, NOT A RESULT. abr-tuning-levers records that a
% +-2% central difference is not trustworthy alone here, and that a lever was
% nearly discarded on one. So this SWEEPS along u and reports what actually
% happens, including the two columns the projection does not control -- maperr
% and oscillatory margin -- either of which can veto the direction.
PROJ = '/Users/neely/Library/CloudStorage/OneDrive-FatherFlanagan''sBoysHome/STN/mod/mccm';
addpath(PROJ); cd(PROJ);
A=load('/Users/neely/mccm_runs/latatlas.mat');
L=load('fit_nch3_surface.mat'); pa=L.R.pa;
if (isfield(pa,'hbmode')), pa=rmfield(pa,'hbmode'); end
H=parfit26('handles'); pv0=H.getpar(pa);

gs=A.ds(:)'; gl=A.dl(:)'; lbl=A.lbl;
live = isfinite(gs) & isfinite(gl) & (abs(gs)>1e-4);   % drop the numerically dead
gs(~live)=0; gl(~live)=0;
fprintf('\n  live parameters: %d of %d  (dropped |dslope|<1e-4, incl. the inert *q set)\n', ...
        sum(live), numel(live));

u = gs - (dot(gs,gl)/max(dot(gl,gl),eps))*gl;    % project out the level_c direction
u(~live)=0;
u = u/max(abs(u));                               % max 1 unit = one 5% step
dsdir = dot(gs,u); dldir = dot(gl,u);
fprintf('  directional derivative:  dslope %+.5f  dlvlc %+.5f  per unit step\n', dsdir, dldir);
fprintf('  (level_c residual is %.1e of the slope movement -- the projection worked)\n', ...
        abs(dldir)/max(abs(dsdir),eps));
need = (0.413 - A.m0.slope)/dsdir;
fprintf('  slope %.4f -> 0.413 needs step s = %+.2f\n', A.m0.slope, need);
[~,ord]=sort(abs(u),'descend');
fprintf('\n  top contributors to the direction:\n');
for k=1:8
    i=ord(k); if (u(i)==0), break; end
    fprintf('    %-8s %+7.3f\n', lbl{i}, u(i));
end

ss = unique(round([0 need/4 need/2 3*need/4 need 1.25*need 1.5*need],3));
ss = sort(ss);
fprintf('\n  %8s %8s %8s %9s %9s %9s\n','step','slope','lvl_c','shoulder','maperr','osc');
res=[];
for i=1:numel(ss)
    s=ss(i); pvk=pv0;
    for j=1:numel(u)
        if (u(j)==0), continue; end
        if (abs(pv0(j))>1e-12), pvk(j)=pv0(j)*(1+0.05*s*u(j));
        else,                   pvk(j)=pv0(j)+0.05*s*u(j); end
    end
    p=H.setpar(pa,pvk,[]);
    try
        m=abr_metric(p,false);
        if (~m.ok), fprintf('  %8.3f  model not ok: %s\n', s, m.msg); continue; end
        Rf=fdm26(struct('pa',p));
        evalc('C=tdm26(''coupeig'',struct(''pa'',p));');
        fprintf('  %8.3f %8.4f %8.3f %9.4f %9.2f %9.1f\n', ...
                s, m.slope, m.level_c, m.shoulder, Rf.maperr, C.maxRe_osc);
        res(end+1,:)=[s m.slope m.level_c m.shoulder Rf.maperr C.maxRe_osc]; %#ok<SAGROW>
    catch e
        fprintf('  %8.3f  FAILED: %s\n', s, e.message(1:min(50,end)));
    end
end
if (size(res,1)>2)
    fprintf('\n  LINEARITY CHECK: predicted slope at s=%.2f is %.4f, measured %.4f\n', ...
            res(end,1), A.m0.slope+dsdir*res(end,1), res(end,2));
    fprintf('  level_c drift over the sweep: %.3f -> %.3f (band [3.5 6.5])\n', ...
            res(1,3), res(end,3));
    good = res(res(:,6)<-40 & res(:,5)<150 & res(:,4)<=res(1,4)+0.02, :);
    if (~isempty(good))
        [~,j]=min(abs(good(:,2)-0.413));
        fprintf('  BEST admissible (osc<-40, maperr<150, shoulder not worse):\n');
        fprintf('    step %.3f -> slope %.4f, lvl_c %.3f, shldr %.4f, maperr %.2f, osc %.1f\n', ...
                good(j,1), good(j,2), good(j,3), good(j,4), good(j,5), good(j,6));
    else
        fprintf('  NO admissible point -- the direction is vetoed by maperr/osc/shoulder.\n');
    end
end
save('/Users/neely/mccm_runs/nulldir.mat','res','u','lbl','dsdir','dldir');
disp('NULLDIR_DONE');
