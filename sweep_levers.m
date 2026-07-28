function sweep_levers
% Lever sweeps, v2. All pa-field changes only (no tdm26 edits), so a concurrent
% parfit26 is unaffected.
%
% CHANGES after v1 was killed at point 2:
%  * PRE-SCREEN the cheap constraints (fdm26 maperr, coupeig stability) BEFORE
%    the expensive abr_metric (16 tbabr runs). v1 wasted a full sweep on a point
%    that was already unstable (osc=+1701) and had maperr 1957.
%  * SAVE INCREMENTALLY after every point (v1 saved only at the end and lost all).
%  * Sweep A step sizes cut ~4x: v1's s=-0.60 destabilized. NOTE v1 FALSIFIED the
%    premise that shifting k1e,m1e together is CF-map-preserving -- CF ~ sqrt(k1/m1)
%    holds for the ISOLATED BM resonance, but the model's CF map also depends on
%    m1 individually via the fluid coupling abmom = bw*dx/m1. So A is now a
%    COST-per-slope probe: how much maperr must we pay per unit of slope movement?

L=load('parfit26_recip.mat'); pa0=L.R.pa;
fprintf('baseline chsz=[%s]  k1e=%.4f m1e=%.4f\n', num2str(pa0.chsz), pa0.k1e, pa0.m1e);
res=[];
res=evalpt(pa0,'baseline',0,res);

fprintf('\n-- A: k1e,m1e co-shift (cost-per-slope probe; NOT map-free) --\n');
for s=[-0.20 -0.10 0.10 0.20]
    p=pa0; p.k1e=pa0.k1e+s; p.m1e=pa0.m1e+s;
    res=evalpt(p,'A_kmshift',s,res);
end

fprintf('\n-- B: place-dependent CA gain gampro=exp(g*(x/L-0.5)) --\n');
n=pa0.n; xn=((1:n)'-1)/(n-1);
for g=[-0.8 -0.4 0.4 0.8]
    p=pa0; p.gampro=exp(g*(xn-0.5));
    res=evalpt(p,'B_gampro',g,res);
end

fprintf('\n-- C: chamber sizes (scale middle/SS chamber; fitted value 0.42) --\n');
for mult=[0.25 0.5 2 4]
    p=pa0; p.chsz=pa0.chsz; p.chsz(2)=pa0.chsz(2)*mult;
    res=evalpt(p,'C_chsz2',mult,res);
end

report(res);
disp('SWEEP_LEVERS_DONE');
end

% -------------------------------------------------------------------------
function res=evalpt(pa,lever,val,res)
k=numel(res)+1;
res(k).lever=lever; res(k).val=val;
res(k).slope=NaN; res(k).level=NaN; res(k).anchor=NaN; res(k).shoulder=NaN;
res(k).maperr=NaN; res(k).osc=NaN; res(k).nsub=NaN; res(k).why='';
fprintf('  %-10s %+6.2f : ',lever,val);
try
    % ---- cheap constraints first -------------------------------------
    Rf=fdm26(struct('pa',pa)); mp=Rf.maperr;
    evalc('S=tdm26(''coupeig'',struct(''pa'',pa));'); osc=S.maxRe_osc;
    res(k).maperr=mp; res(k).osc=osc;
    if (~isfinite(mp) || mp>500 || ~isfinite(osc) || osc>0)
        res(k).why=sprintf('screened (maperr=%.0f osc=%+.0f)',mp,osc);
        fprintf('%s -- skipped tbabr\n',res(k).why);
        save('sweep_levers.mat','res'); return
    end
    % ---- expensive ABR surface ---------------------------------------
    m=abr_metric(pa,false);
    if (~m.ok), res(k).why=['abr_metric: ' m.msg]; fprintf('FAILED %s\n',m.msg);
        save('sweep_levers.mat','res'); return; end
    fi=find(abs(m.f-1)<0.01,1); li=find(abs(m.slv-20)<0.1,1);
    res(k).slope=m.slope; res(k).level=100*(m.level_c^(1/100)-1);
    res(k).anchor=m.lat(fi,li); res(k).shoulder=m.shoulder; res(k).nsub=m.n_sub;
    [~,D]=abr_surface_obj(m);            % full-surface view too
    res(k).b=D.b; res(k).dsurf=D.d; res(k).shape=D.resid;
    fprintf('slope=%.3f b=%.2f shape=%.3f maperr=%.0f shoul=%.3f osc=%+.0f\n', ...
            m.slope,D.b,D.resid,mp,m.shoulder,osc);
catch e
    res(k).why=['ERROR ' e.message]; fprintf('ERROR %s\n',e.message);
end
save('sweep_levers.mat','res');
end

% -------------------------------------------------------------------------
function report(res)
fprintf('\n%-11s %6s %7s %7s %7s %8s %8s %8s %7s\n', ...
        'lever','val','slope','b','shape','shoulder','maperr','osc','anchor');
fprintf('%s\n',repmat('-',1,80));
for k=1:numel(res)
    r=res(k);
    if (~isempty(r.why) && isnan(r.slope))
        fprintf('%-11s %6.2f   %s\n', r.lever, r.val, r.why); continue
    end
    b=NaN; sh=NaN;
    if (isfield(r,'b')), b=r.b; end
    if (isfield(r,'shape')), sh=r.shape; end
    fprintf('%-11s %6.2f %7.3f %7.2f %7.3f %8.3f %8.0f %8.0f %7.2f\n', ...
        r.lever,r.val,r.slope,b,sh,r.shoulder,r.maperr,r.osc,r.anchor);
end
fprintf(['\ntargets: slope 0.39-0.41 | b 11.99-13.27 (baseline 11.5 = TOO FAST) | ' ...
         'shape->0 | shoulder->0 | maperr<=185 | osc<0\n']);
end
