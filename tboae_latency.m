function [lat, emag, out] = tboae_latency(pr, ps, fr, dtms)
%TBOAE_LATENCY  Latency (energy-weighted group delay) of a tone-burst OAE.
%   [lat,emag,out] = tboae_latency(pr, ps, fr, dtms)
%
%   pr   : eardrum-pressure trace WITH roughness present (the physical cochlea)
%   ps   : eardrum-pressure trace of the SMOOTH reference (rough_amp = 0)
%   fr   : stimulus (tone-burst) frequency, kHz
%   dtms : sample interval, ms
%
%   The emission is  e = pr - ps.  Differencing cancels the identical forward
%   stimulus drive and the passive middle-ear response, leaving ONLY the
%   roughness-induced coherent-reflection emission -- so no stimulus-artifact
%   gating is required (unlike a raw ear-canal recording).  Cf. fdm12a Fig 7
%   (negligible reflection when smooth; significant reflection when rough).
%
%   lat  : energy-weighted centroid time of the on-frequency emission envelope
%          (ms).  This is the tone-burst OAE latency to compare against the
%          2013 power-law fit and the ~2x-forward round-trip ratio.
%   emag : rms level of the band-limited emission (dB, arbitrary common ref).
%   out  : struct .e (band emission) .env (envelope) .t (ms) .gd (spectral
%          energy-weighted group delay, ms, cross-check) .snr (env peak/rms).
%
%   Called with NO arguments, runs a self-test on a synthetic delayed burst.

if (nargin == 0), selftest; return; end

pr = pr(:); ps = ps(:); n = numel(pr);
e  = pr - ps; e = e - mean(e);                 % emission (stimulus cancels)

fs = 1000/dtms;                                % sample rate (Hz)
f0 = fr*1000;                                  % stimulus frequency (Hz)

% ---- signed bin frequencies and a half-octave analysis band around f0 -------
k  = (0:n-1)';
fk = k*(fs/n); fk(fk > fs/2) = fk(fk > fs/2) - fs;   % signed freq (Hz)
inband = abs(fk) >= f0/sqrt(2) & abs(fk) <= f0*sqrt(2);

% ---- band-limited analytic signal (manual Hilbert; no toolbox needed) -------
E  = fft(e);
Za = zeros(n,1);
pos = fk > 0 & inband;                         % keep positive in-band, x2
Za(pos) = 2*E(pos);
z   = ifft(Za);                                % complex analytic band signal
env = abs(z);                                  % emission envelope
eb  = real(z);                                 % band-limited real emission

% ---- latency = envelope-peak arrival, refined by a peak-anchored energy-
%      weighted centroid (a robust "group delay").  A WHOLE-record centroid is
%      deliberately NOT used: when the emission is weak (below the model's
%      numerical floor) it collapses to the record center (~T/2) for every
%      condition -- an artifact.  Peak-anchoring localizes the real arrival.
t    = (0:n-1)'*dtms;                          % ms
tlo  = 1; thi = min(35, t(end)-1); wcw = 3;    % analysis window (ms) & centroid half-width
emag = 20*log10(sqrt(mean(eb.^2)) + eps);      % band emission level (dB)
win  = t >= tlo & t <= thi;
if (~any(win) || ~all(isfinite(env)) || max(env(win)) <= 0)
    lat = NaN; tpk = NaN;
else
    ew = env; ew(~win) = 0;
    [~,ip] = max(ew); tpk = t(ip);             % envelope-peak arrival
    wc = t >= tpk-wcw & t <= tpk+wcw;          % local, peak-anchored window
    lat = sum(t(wc).*env(wc).^2)/sum(env(wc).^2);  % energy-weighted group delay (ms)
end

if (nargout > 2)
    out.e = eb; out.env = env; out.t = t; out.tpk = tpk; out.emag = emag;
    out.reclen = t(end);
    out.snr = max(env)/(sqrt(mean(env.^2)) + eps);
    % spectral energy-weighted group delay (one-sided), independent cross-check
    np   = floor(n/2);
    Ep   = E(1:np+1); fp = (0:np)'*(fs/n);
    phi  = unwrap(angle(Ep));
    gd   = -gradient(phi, 2*pi*fp)*1000;       % ms
    bp   = fp >= f0/sqrt(2) & fp <= f0*sqrt(2);
    Wp   = abs(Ep).^2;
    if (any(bp) && sum(Wp(bp)) > 0)
        out.gd = sum(gd(bp).*Wp(bp))/sum(Wp(bp));
    else
        out.gd = NaN;
    end
end
end % tboae_latency

% -------------------------------------------------------------------------
function selftest
fs = 100000; dtms = 1000/fs; T = 0.050; n = round(T*fs);
t  = (0:n-1)'*dtms;                            % ms
fr = 2; f0 = fr*1000;                          % 2 kHz burst

T0 = 8; sig = 2.0;                             % true emission delay & width (ms)
carrier = sin(2*pi*f0*t/1000);
stim = 5*carrier .* (t < 2);                   % common stimulus part (early)
emis = carrier .* exp(-((t - T0)/sig).^2);     % delayed reflection emission
pr = stim + 0.01*emis;                         % rough: stimulus + emission
ps = stim;                                     % smooth: stimulus only

[lat, emag, out] = tboae_latency(pr, ps, fr, dtms);
fprintf(['tboae_latency self-test: true delay=%.2f ms  centroid=%.2f ms  ' ...
         'spectral-gd=%.2f ms  emag=%.1f dB  snr=%.1f\n'], ...
         T0, lat, out.gd, emag, out.snr);
ok1 = abs(lat - T0)   < 0.5;
ok2 = abs(out.gd - T0) < 0.8;
fprintf('  centroid within 0.5 ms: %d   spectral-gd within 0.8 ms: %d\n', ok1, ok2);
if (ok1 && ok2), fprintf('  SELF-TEST PASS\n'); else, fprintf('  SELF-TEST FAIL\n'); end
end
