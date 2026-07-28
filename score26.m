function S = score26(pa, mode, verbose)
%SCORE26  Unified scorecard for any mccm configuration, chambers 1-4.
%
%   S = score26(pa)            full scorecard (TBABR + TBOAE + fast block)
%   S = score26(pa,'fast')     click-based block only (no TBABR/TBOAE)
%
%   PURPOSE. Persistent incremental improvement requires a COMMON YARDSTICK:
%   results scored by differently-implemented metrics cannot be compared, so
%   gains cannot be banked. Metric REIMPLEMENTATION has been this project's
%   dominant failure mode (the peak-switching latency detector; mass_cost.m
%   wrong in both columns via sparse-vs-dense place sampling; tiptail_metric's
%   0.6-16 kHz window silently rejecting every m=4 model). score26 therefore
%   REUSES the already-validated routines and reimplements nothing:
%       abr_metric-style tbabr run -> tdm26('tbabr',...)   (latency surface)
%       abr_surface_obj                                    (surface fit b,c,d)
%       tboae_latency via pa.oae                           (emission)
%       fdm26                                              (maperr / tuning)
%       tdm26('coupeig')                                   (stability)
%   The click block below is the local_tip measurement validated in
%   ohcgain_sweep.m against a known control (native m=4 gain = +2.40 dB).
%
%   TWO CHAMBER-COUNT-DEPENDENT SUBTLETIES, both previously cost real time:
%   (1) THE AMPLITUDE KNOB DIFFERS. For m<4, pa.gam scales the active term and
%       the model is passive at gam=0. For m>=4, gam scales ONLY the k4/r4
%       sub-term -- the gh*k3 active coupling SURVIVES at gam=0 -- so the whole
%       OHC force is scaled by pa.ohcgain instead. The passive reference is
%       chosen accordingly; using gam for m=4 measures the wrong thing.
%   (2) NO HARD-CODED CF WINDOW. The BF map is scored over whatever band the
%       model actually occupies (m=4 spans ~0.3-2 kHz natively, m=2 ~0.3-16),
%       so a 2-chamber-shaped window must never gate a 4-chamber model.
%
%   Latency is measured with hbmode='bm' (BM-peak) -- the corrected metric. The
%   default mixed drive peak-switches between two physical peaks as level rises
%   and manufactured the discredited 0.413 result.
%
%   Components are reported SEPARATELY and deliberately NOT combined into one
%   weighted scalar: how TBABR and TBOAE should trade off is a scientific
%   choice, not a default.

if (nargin<2 || isempty(mode)), mode='full'; end
if (nargin<3), verbose=true; end

S = struct('mode',mode,'nch',numel(pa.chsz),'m',pa.m,'when',datestr(now)); %#ok<TNOW1,DATST>
% Place grid as FRACTIONS of pa.n, not fixed indices. The original
% 1391:-10:11 was valid only for n=1401; at any other grid size those
% indices fall outside the array. Fractions keep the SAME physical places
% under a change of n, which is what makes a convergence test meaningful.
ISVFRAC = (1391:-10:11)/1401;
ISV = unique(max(1, min(pa.n, round(ISVFRAC*pa.n))), 'stable');
XLO = 0.05; XHI = 0.85; FP = [1 2 4 8];

% ================= FAST BLOCK: click-derived map / tip / amplifier ==========
pa1 = pa; pa0 = pa;
if (pa.m >= 4)
    % ACTIVE level = the config's OWN ohcgain (default 1 if unset). Forcing 1
    % here made amp_gain blind to a swept ohcgain: the c4_reasonable sweep
    % returned +2.46 dB identically for ohcgain 1, 2 and 4 because every row
    % was measured at 1. Passive reference stays 0.
    og = 1; if (isfield(pa,'ohcgain')), og = pa.ohcgain; end
    pa1.ohcgain = og;  pa0.ohcgain = 0;
else
    pa1.gam = pa.gam; pa0.gam = 0;      % gam is the clean knob below m=4
end
t1 = local_tip(pa1, ISV, XLO, XHI, FP);
t0 = local_tip(pa0, ISV, XLO, XHI, FP);
S.click_ok   = strcmp(t1.fin,'yes');
S.click_note = t1.fin;
S.bf_range   = t1.rng;        % tonotopic span, octaves (m=2 native ~5.9)
S.bf_mono    = t1.mono;       % strictly increasing BF apex->base?
S.bf_lo      = t1.bflo; S.bf_hi = t1.bfhi; S.bf_fold = t1.fold;
S.contrast   = t1.chi;        % tip-tail contrast at the most basal valid place
S.amp_gain   = t1.pk(2) - t0.pk(2);   % active-vs-passive gain at 2 kHz, dB (BM)
S.amp_d2     = t1.pk2(2) - t0.pk2(2); % same for the shear coordinate
S.amp_d3     = t1.pk3(2) - t0.pk3(2); % same for OC height (m>=4 only)
S.amp_max    = max([S.amp_gain S.amp_d2 S.amp_d3]);  % largest of the three
S.xbest      = t1.xb;  S.lat_click = t1.lt;

% ---- tuning / MAP (frequency domain; independent of the latency detector) --
S.maperr = NaN;
try
    Rf = fdm26(struct('pa',pa));
    if (isstruct(Rf) && isfield(Rf,'maperr')), S.maperr = Rf.maperr; end
catch e
    S.maperr_err = e.message;
end

% ---- stability: report BOTH. maxRe_osc filters to oscillatory modes >1 kHz
%      and is BLIND to static (f=0) divergence -- it once read -486 on a model
%      whose time march blew up to 3e296. maxRe(all) is the real arbiter.
S.maxRe = NaN; S.maxRe_osc = NaN;
try
    evalc('E = tdm26(''coupeig'', struct(''pa'',pa));');
    S.maxRe = E.maxRe; S.maxRe_osc = E.maxRe_osc;
catch e
    S.coupeig_err = e.message;
end

if (strcmpi(mode,'fast')), if (verbose), print_card(S); end, return; end

% ================= FULL BLOCK: TBABR latency surface + TBOAE ================
% One tbabr run yields both: pa.oae makes the tone-burst OAE a secondary output
% of the SAME stimuli (so ABR and OAE latencies are directly comparable), at
% roughly 2x cost because each condition also runs a smooth (rough_amp=0)
% reference for the differencing.
pr = pa;
pr.hbmode = 'bm';                       % BM-peak latency (corrected metric)
pr.oae    = 1;
if (~isfield(pr,'rough_amp') || pr.rough_amp < 1e-2)
    pr.rough_amp = 3e-2;                % below ~1e-2 the emission sits under
end                                     % the numerical floor -> flat latency
S.rough_amp = pr.rough_amp;
S.tbabr_ok = false;
try
    evalc('T = tdm26(''tbabr'', pr, 0, 0);');
    S.tbabr_ok = true;
catch e
    S.tbabr_err = e.message;
    if (verbose), print_card(S); end
    return
end
S.f = T.f(:).'; S.slv = T.slv(:).'; S.lat = T.lat; S.abr = T.abr;
S.n_sub = nnz(~isfinite(T.lat));
if (isfield(T,'sho')), S.sho = T.sho; S.shoulder = mean(T.sho(isfinite(T.sho))); end

% ---- TBABR: full 16-cell surface fit to tau = b*c^(-i)*f^(-d) --------------
try
    m = struct('f',S.f,'slv',S.slv,'lat',S.lat);
    [J,D] = abr_surface_obj(m);
    S.surf_J     = J;        % objective (shape + hinge band losses)
    S.surf_resid = D.resid;  % rms log-residual = SHAPE error, band-independent
    S.b = D.b; S.c = D.c; S.d = D.d;
    S.hinge_d = D.hd; S.hinge_c = D.hc; S.hinge_b = D.hb;
    S.surf_n  = D.n;
    S.level_pct = 100*(D.c^(1/100) - 1);   % %/dB, the commensurate level metric
catch e
    S.surf_err = e.message;
end

% ---- TBOAE: emission latency / magnitude, and the round-trip ratio ---------
if (isfield(T,'oae')), S.oae = T.oae; end
if (isfield(T,'oam')), S.oam = T.oam; end
if (isfield(S,'oae') && ~isempty(S.oae))
    o = S.oae; l = S.lat;
    ok = isfinite(o) & isfinite(l) & l > 0;
    S.oae_n = nnz(ok);
    if (S.oae_n > 0)
        S.oae_lat  = mean(o(ok));
        S.oae_ratio = mean(o(ok) ./ l(ok));   % expect ~2 (round trip) if the
    end                                       % emission is coherent reflection
end
if (isfield(S,'oam')), S.oae_mag = mean(S.oam(isfinite(S.oam))); end

if (verbose), print_card(S); end
end % score26

% ---------------------------------------------------------------------------
function print_card(S)
fprintf('\n=== score26  m=%d (%d chambers)  [%s] ===\n', S.m, S.nch, S.mode);
fprintf('  MAP/tuning   maperr    %s\n', num(S,'maperr','%8.1f'));
fprintf('  Stability    maxRe     %s   maxRe_osc %s\n', ...
        num(S,'maxRe','%+8.1f'), num(S,'maxRe_osc','%+8.1f'));
fprintf('  CF map       range     %s oct   [%s .. %s kHz]  mono=%s\n', ...
        num(S,'bf_range','%5.2f'), num(S,'bf_lo','%.2f'), num(S,'bf_hi','%.2f'), ...
        tern(isfield(S,'bf_mono'), S_get(S,'bf_mono'), '?'));
fprintf('               fold %s oct (largest downward BF step; ~0 = clean)\n', num(S,'bf_fold','%5.2f'));
fprintf('  Tip          contrast  %s dB   amp(active-passive) %s dB\n', ...
        num(S,'contrast','%6.1f'), num(S,'amp_gain','%+6.2f'));
fprintf('  Amp by DOF   d1(BM) %s   d2(shear) %s   d3(OC ht) %s\n', ...
        num(S,'amp_gain','%+6.2f'), num(S,'amp_d2','%+6.2f'), num(S,'amp_d3','%+6.2f'));
if (isfield(S,'surf_resid'))
    fprintf('  TBABR surf   J %s  shape-rms %s   n=%s\n', ...
            num(S,'surf_J','%6.3f'), num(S,'surf_resid','%6.3f'), num(S,'surf_n','%d'));
    fprintf('  TBABR law    b %s ms (tgt 11.99-13.27)  c %s (tgt 5.0-5.34)  d %s (tgt .39-.41)\n', ...
            num(S,'b','%6.2f'), num(S,'c','%5.2f'), num(S,'d','%5.3f'));
    fprintf('               level %s %%/dB (tgt 1.62)   hinges d/c/b %s %s %s\n', ...
            num(S,'level_pct','%5.2f'), num(S,'hinge_d','%.3f'), ...
            num(S,'hinge_c','%.3f'), num(S,'hinge_b','%.3f'));
end
if (isfield(S,'oae_lat'))
    fprintf('  TBOAE        lat %s ms  ratio-to-ABR %s (expect ~2)  mag %s dB  n=%s\n', ...
            num(S,'oae_lat','%6.2f'), num(S,'oae_ratio','%5.2f'), ...
            num(S,'oae_mag','%7.1f'), num(S,'oae_n','%d'));
end
if (isfield(S,'n_sub')), fprintf('  sub-threshold cells: %d\n', S.n_sub); end
if (isfield(S,'tbabr_err')), fprintf('  TBABR FAILED: %s\n', S.tbabr_err); end
if (~S.click_ok), fprintf('  CLICK PROBLEM: %s\n', S.click_note); end
end

function s = num(S,f,fmt)
if (isfield(S,f) && ~isempty(S.(f)) && isnumeric(S.(f)) && isfinite(S.(f)(1)))
    s = sprintf(fmt, S.(f)(1));
else
    s = '     --';
end
end
function v = S_get(S,f), v = S.(f); end
function s = tern(c,a,b), if c, s=a; else, s=b; end, end

% ---------------------------------------------------------------------------
function t = local_tip(pa, ISV, XLO, XHI, fp)
% VERBATIM the measurement validated in ohcgain_sweep.m (native m=4 gain
% +2.40 dB control), extended to also return the BF-map range. Interior-only
% scoring: the dense grid otherwise reaches the helicotrema/stapes and the peak
% search lands on boundary artifacts.
t.fin='yes'; t.chi=NaN; t.deg=1; t.rng=NaN; t.mono='?'; t.bflo=NaN; t.bfhi=NaN; t.fold=NaN;
t.xb=nan(1,numel(fp)); t.lt=t.xb; t.pk=t.xb;
% PEAKS IN ALL THREE PARTITION COORDINATES. Amplification need not appear in
% the BM: species with no basilar membrane (lizards, birds, amphibians) show
% sharp tuning and active emissions, so the amplifier is a property of the hair
% cell and tectorial micromechanics, and the BM's distinct role is the
% tonotopic travelling wave. Measuring gain on d1 alone can therefore MISS a
% working amplifier that acts locally on the shear (d2) or OC height (d3) and
% couples only weakly to the BM. m=1/m=2 have no separate micromechanical
% coordinates, so all their gain necessarily appears in d1 by construction --
% which means a d1-only comparison across chamber counts is not comparing the
% same physical quantity.
t.pk2=t.xb; t.pk3=t.xb;
pa.isv=ISV;
try, evalc('S=tdm26(0,pa,0,0);'); catch, t.fin='THREW'; return; end
if (any(~isfinite(S.d1(:)))), t.fin='DIVERG'; return; end
nf=numel(S.f); f=S.f(:); P=fft(S.ped); P=P(1:nf); np=size(S.d1,2);
xp=pa.isv(:)'/pa.n; inr=xp>=XLO & xp<=XHI;
H=zeros(nf,np); for i=1:np, D=fft(S.d1(:,i)); H(:,i)=D(1:nf)./max(abs(P),eps); end
Am=abs(H); co=nan(1,np); bf=nan(1,np);
% same transfer for d2 and d3 (d3 is all zeros below m=4, handled by the guard)
Am2=zeros(nf,np); Am3=zeros(nf,np);
for i=1:np
    D2=fft(S.d2(:,i)); Am2(:,i)=abs(D2(1:nf)./max(abs(P),eps));
    if (isfield(S,'d3') && ~isempty(S.d3) && any(S.d3(:)~=0))
        D3=fft(S.d3(:,i)); Am3(:,i)=abs(D3(1:nf)./max(abs(P),eps));
    end
end
% ---- CONTINUITY-CONSTRAINED BF TRACKER (2026-07-27) ----------------------
% REPLACES a per-place global argmax, which had two failure modes.
%
% (1) EDGE. Am is f-WEIGHTED, so where the response above CF is negligible (the
% apical places) the weighted curve rises monotonically to the band edge and the
% argmax lands ON it. tdm26's unbanded find_bf returned the top bin outright --
% two apical places both reporting 24.98 kHz at -39 dB. An edge maximum means
% the curve never turned over, so no BF exists there and NaN is correct.
%
% (2) MODE COMPETITION, which the edge bug was masking. fold_probe measured the
% residual 0.84 oct fold in m=3b directly: at x/L 0.101 the weighted curve peaks
% at 10.449 kHz with a second peak at 5.200 kHz just 0.29 dB down; one place
% further basal, at x/L 0.094, the SAME pair reads 5.835 kHz and 10.791 kHz
% separated by 0.01 dB. The peaks swap rank, and a global argmax reports a
% 0.84 oct cliff where the mechanics are perfectly smooth. That is the BM mode
% and the shear mode crossing -- their separation is ~1.0 oct at that place,
% matching the step. m=3b's map was never defective; the scalar summary was.
%
% THE TRACKER seeds at the place whose dominant peak is CLEAREST (largest margin
% over its runner-up, so the seed is the one place least likely to be
% mis-assigned) and walks outward in both directions, taking at each step the
% local maximum nearest IN LOG FREQUENCY to the previously accepted BF. A step
% larger than MAXJUMP is refused and left NaN rather than followed, so one bad
% place cannot derail the rest of the track.
%
% CONSEQUENCE t.rng, t.fold and t.mono now describe the MAP. Previously they
% described the mode ORDERING, which is why FOLD fired on every model in the
% project including m=3b, the accepted scaffolding.
% A CONTINUITY-CONSTRAINED TRACKER WAS TRIED HERE AND REVERTED. Do not re-add
% it without reading this. The reasoning that motivated it was wrong.
%
% Continuity and CF are DIFFERENT QUANTITIES. At a mode crossing the mode
% trajectory continues smoothly while the DOMINANT peak switches branches, so a
% continuity tracker follows the weaker mode past the crossing BY DESIGN. CF is
% defined by the largest response, so the argmax jump is not a detector failure.
%
% Measured (tracker_diag.m, 2026-07-27, m=4 new default): the tracker left the
% dominant mode at 15 of 106 valid places, x/L 0.051 to 0.151, in ONE CONTIGUOUS
% RUN beginning at the crossing, with the deficit growing monotonically to
% -16.21 dB and the tracked peak falling to RANK 11. It also drove contrast
% NEGATIVE (-13.6 vs +6.4) and doubled bf_hi (7.40 -> 14.77 kHz), both direct
% consequences of reporting a peak that is not the largest.
%
% WHAT THE FOLD ACTUALLY MEASURES: mode degeneracy, not map quality.
% fold_probe.m found the m=3b fold at x/L 0.101/0.094 where the BM mode and the
% shear mode cross with the runner-up just 0.29 dB and then 0.01 dB down. At
% that margin the CF is genuinely knife-edge. A residual fold is a REAL property
% of the model's mode structure and should be read as such -- it is not a defect
% to engineer away, and m=3b's map was never broken.
%
% The EDGE guard below is a genuine bug fix and stays.
for i=find(inr)
    a=Am(:,i).*f; bnd=(f>0.15&f<18); a(~bnd)=0; [q,ip]=max(a); if(q<=0),continue;end
    % APICAL EDGE GUARD (2026-07-27), same fault as tdm26's find_bf. Am is
    % f-WEIGHTED, so at apical places, where the response above CF is
    % negligible, the weighted curve rises monotonically to the edge of the band
    % and the argmax lands ON that edge. tdm26's unbanded version returned the
    % top bin outright (two apical places both reporting 24.98 kHz at -39 dB);
    % banded here the same failure pins bf to the 18 kHz limit instead, which is
    % harder to spot. An edge maximum means the curve never turned over, so no
    % BF exists at this place. Leaving bf(i) NaN is correct -- callers filter on
    % isfinite, and a wrong BF corrupts rng, fold and mono for the whole array.
    % This discards ~25% of places, all apical, so rng UNDERSTATES the true span.
    ib=find(bnd);
    if (isempty(ib) || ip<=ib(1) || ip>=ib(end)), continue; end
    bf(i)=f(ip);
    [~,it]=min(abs(f-f(ip)/4)); co(i)=20*log10(q/max(a(it),eps));
end
ii=find(isfinite(co)); if(~isempty(ii)), t.chi=co(ii(end)); end
v=bf(isfinite(bf)&bf>0);
if (numel(v)>=2)
    t.rng=log2(max(v)/min(v)); t.bflo=min(v); t.bfhi=max(v);
    % FOLD as a MAGNITUDE, not a binary flag. On 139 dense places even the
    % m=2 native reference (a good, textbook map) shows sub-resolution BF
    % wiggle, so all(diff>0) trips on the gold standard and cannot
    % discriminate. t.fold = largest single DOWNWARD BF step, in octaves,
    % apex->base: ~0 for a clean map, ~0.4 for the m=4 mid-array reversal.
    dd = diff(log2(v));
    t.fold = max([0, -min(dd)]);
    t.mono = tern(t.fold < 0.10, 'ok', 'FOLD');
end
w=2*pi*f*1000;
for k=1:numel(fp)
    [~,jf]=min(abs(f-fp(k))); e2=Am(jf,:); e2(~isfinite(e2))=0; e2(~inr)=0;
    [~,jb]=max(e2); if(e2(jb)<=0),continue;end
    t.xb(k)=xp(jb); t.pk(k)=20*log10(e2(jb));
    q2=Am2(jf,:); q2(~isfinite(q2))=0; q2(~inr)=0; m2=max(q2);
    q3=Am3(jf,:); q3(~isfinite(q3))=0; q3(~inr)=0; m3=max(q3);
    if (m2>0), t.pk2(k)=20*log10(m2); end
    if (m3>0), t.pk3(k)=20*log10(m3); end
    p2=unwrap(angle(H(:,jb))); g2=-gradient(p2,w);
    t.lt(k)=median(g2(max(1,jf-3):min(nf,jf+3)))*1000;
end
t.deg=(std(t.xb,0,'omitnan')<1e-3)||(median(abs(t.lt),'omitnan')<0.30);
end
