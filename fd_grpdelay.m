function R = fd_grpdelay(pa, flst)
% FD_GRPDELAY  Frequency-domain cochlear group delay (forward-latency analog).
%   Thin wrapper around fdm26's group-delay analysis (reuses the internal
%   fdmod23 solver). Returns R.f (CF, kHz), R.tau (ms), R.d (log-log slope over
%   0.5-4 kHz; ABR target 0.413), R.tipgain (max active tuning gain, dB -- a
%   near-criticality proxy).
%
%   R = FD_GRPDELAY(PA[, FLST])  PA from modpar26 (fields may be modified to
%   probe dispersion params k1e/m1e/xtap/xtex); FLST in Hz (default fine grid).
pr.pa = pa;
if (nargin>1 && ~isempty(flst)), pr.flst = flst; end
R = fdm26(pr);
end
