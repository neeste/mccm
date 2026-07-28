function [J,D] = abr_surface_obj(m, opts)
%ABR_SURFACE_OBJ  Full 16-point latency-surface objective with band losses.
%
%   Replaces the summary-statistic objective (which compresses all 16 measured
%   latencies into just slope + level_c, discarding most of the information and
%   flattening the search landscape) with a weighted fit of the FULL surface to
%
%        tau = b * c^(-i) * f^(-d),      i = dB/100
%
%   plus HINGE (band) losses that are exactly ZERO inside the 2013 replication's
%   published uncertainty -- so the fit is not driven to spurious precision and
%   does not trade real quality for decimal places on a point target.
%
%   [J,D] = abr_surface_obj(m,opts)
%     m    : abr_metric output (uses .f, .slv, .lat)
%     opts : .dband [0.39 0.41]  .cband [5.0 5.34]  .bband [11.99 13.27]
%            .wshape .wd .wc .wb
%   D     : diagnostics -- fitted b,c,d, rms log-residual (SHAPE error), mask.
%
%   The rms log-residual is information the old objective never saw: it measures
%   whether the model's latency surface is power-law SHAPED at all, independent
%   of whether the exponent is right.

if (nargin<2), opts=struct; end
g=@(fn,dv) getdef(opts,fn,dv);
dband=g('dband',[0.39 0.41]); cband=g('cband',[5.0 5.34]);
bband=g('bband',[11.99 13.27]);
wshape=g('wshape',1.0); wd=g('wd',1.0); wc=g('wc',0.5); wb=g('wb',0.25);

f=m.f(:); lv=m.slv(:)'; T=m.lat;
[nf,nl]=size(T);
F=repmat(f,1,nl); I=repmat(lv/100,nf,1);

% ---- mask -------------------------------------------------------------
% The largest-peak detector returns a late (~18-20 ms) value when it finds no
% real onset; those are detector FAILURES, not latencies, and must not be
% fitted. Flag any cell exceeding 3x the nominal 1988-law expectation.
nom = 12.9 .* 5.^(-I) .* F.^(-0.413);
ok  = isfinite(T) & T>0 & T < 3*nom;
D.mask=ok; D.n=sum(ok(:)); D.nom=nom;
if (D.n < 6)
    J=1e6; D.b=NaN; D.c=NaN; D.d=NaN; D.resid=NaN; D.hd=NaN; D.hc=NaN; D.hb=NaN;
    return
end

% ---- weighted LS: log tau = log b - i*log c - d*log f  (linear in params)
y = log(T(ok));
A = [ones(D.n,1), -I(ok), -log(F(ok))];
th = A\y;
D.b = exp(th(1)); D.c = exp(th(2)); D.d = th(3);
r = y - A*th;
D.r = r; D.resid = sqrt(mean(r.^2));      % rms log-residual = SHAPE error

% ---- hinge (band) losses ----------------------------------------------
hin = @(v,b) max(0,b(1)-v) + max(0,v-b(2));
D.hd = hin(D.d,dband); D.hc = hin(D.c,cband); D.hb = hin(D.b,bband);

J = wshape*D.resid + wd*D.hd + wc*D.hc + wb*D.hb;
end

function v=getdef(o,fn,dv)
if (isfield(o,fn)), v=o.(fn); else, v=dv; end
end
