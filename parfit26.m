function R = parfit26(nch, opts)
% PARFIT26  Routine joint fit of the multi-chamber model to MAP (tuning/threshold)
%   AND forward latency (ABR slope), on the current modpar26 / modpar26c3 model.
%
%   Reconciled successor to parfit24. It carries NO duplicated frequency-domain
%   engine: the tuning/MAP error is delegated to fdm26, the forward-latency slope
%   to abr_metric (the time-domain tbabr metric — reliable for the 3-chamber,
%   where the fd group-delay slope is not), and the sub-criticality check to
%   tdm26's coupeig. It fits the impedance parameters AND the chamber sizes (the
%   structural forward-latency lever) together.
%
%     J = wslope*|tdm_slope - 0.413|                       (forward latency, reliable)
%       + wcrit *max(0, coupeig_osc_maxRe + 40)            (sub-critical)
%       + wmap  *max(0, fd_maperr - maptol)                (MAP/tuning fit)
%       + 0.1   *level_c outside [lc(1) lc(2)]             (level dependence)
%
%   R = PARFIT26(NCH, OPTS).  OPTS fields (all optional):
%     .fitidx  params to fit, getpar order 1..30 = impedance, 31..(30+#chsz) =
%              chamber sizes  (default: CF-map + tip + chsz set that closed the fit)
%     .maxfe   fminsearch budget            (default 150)
%     .wslope .wmap .wcrit                  (default 1, 0.001, 0.005)
%     .wstat   weight on the STATIC-divergence guard, max(0, maxRe - statol)
%              (default: wcrit, i.e. unchanged behaviour). Set it independently
%              when the static guard must bite without inflating the osc term --
%              see the note at the J assembly for why one shared weight cannot
%              size both.
%     .chszderive  index of the chsz element DERIVED from the others so that
%              sum(chsz)==2 identically ([] = off, chsz free as before). Set 3
%              to derive SV. Enforced inside setpar_l, so NO search step can
%              violate it -- unlike a post-hoc check, and unlike renormalisation,
%              which was rejected because it makes added chambers resize existing
%              ones. The derived index is auto-removed from fitidx.
%     .maptol  MAP acceptance level         (default 167)
%     .lc      level_c band [lo hi]         (default [3.5 6.5])
%     .warm    starting pa or a .mat filename with R.pa   (default: refit_c3_map.mat
%              if present, else modpar26(nch) baseline)
%     .out     .mat to save R into          (default parfit26_nch<n>.mat)

if (nargin<1 || isempty(nch)), nch=3; end
if (nargin<2), opts=struct(); end
gv=@(f,d) subsref_default(opts,f,d);
fitidx = gv('fitidx', [1 2 3 4 5 7 8 9 10 11 13 20 21 31 32 33]);
% CHAMBER-COUNT AWARE (2026-07-29). pv is [30 impedance, numel(chsz) chamber
% sizes], so indices 31+ address chsz. The default list reaches 33, which is
% valid only for nch>=3: modpar26(1) and modpar26(2) ship chsz=[1 1], giving 32
% entries, so index 33 threw "Index exceeds the number of array elements" and
% nch=1/2 could not be fitted through this path AT ALL with default settings.
% Clip rather than error -- a caller asking for chamber sizes that do not exist
% wants the ones that do.
nc_ = numel(modpar26(nch).chsz);
bad_ = fitidx > (30 + nc_);
if (any(bad_))
    fprintf('  fitidx: dropping %s (nch=%d has only %d chsz entries)\n', ...
            mat2str(fitidx(bad_)), nch, nc_);
    fitidx = fitidx(~bad_);
end
maxfe  = gv('maxfe', 150);
wslope = gv('wslope', 1);  wmap = gv('wmap', 0.001);  wcrit = gv('wcrit', 0.005);
% AMPLIFIER GUARD (2026-07-29). Without it the objective REWARDS deleting the
% cochlear amplifier: m=4 nested-no-vent scores the project's best maperr
% (353.5) with a 4.14-octave map and maxRe 0.0 -- and amp_gain = -0.01 dB. Every
% other term prefers it, so an optimizer will find it reliably.
% wgain*max(0, gainmin - amp_gain) penalises SHORTFALL only, so it cannot reward
% runaway gain (the full 2010 stack hit +1868 dB).
% COST: amp_gain needs a real active-vs-passive comparison, so this calls
% score26('fast') and roughly DOUBLES the per-evaluation cost. wgain=0 disables
% it, but read the next note first.
% fdm26's tipgain was TESTED AS A FREE PROXY AND REJECTED -- it is already
% computed every evaluation, but it ANTI-CORRELATES with gain where it matters:
%     m=1        tipgain 38.56   amp_gain +39.11
%     m=3b       tipgain 50.35   amp_gain +103.02
%     m=4        tipgain 37.87   amp_gain +56.59
%     m=4 no-amp tipgain 62.70   amp_gain -0.01   <-- HIGHEST tipgain of all four
% It tracks tuning SHARPNESS, not active gain, and a passive model can be sharp
% (that config also had the highest tip-tail contrast, 16.07 vs 6.41). Using it
% would have made the objective PREFER the amplifier-free model.
wgain   = gv('wgain', 0.01);  % weight on amplifier shortfall
gainmin = gv('gainmin', 40);  % dB; SN: "OHC forces amplify BM by 40 dB or more"
wlevel = gv('wlevel', 0);    % weight on |level %/dB - 1.62|, %/dB=100*(level_c^(1/100)-1); 0=band only
wanchor= gv('wanchor', 0);   % weight on the ABSOLUTE-latency anchor |WNR(anchor)-target(anchor)| in ms
anchor = gv('anchor', [1 20]);% [freq_kHz level_dB] at which to anchor the WNR latency (target = 1988 formula)
wshoulder = gv('wshoulder', 0); % weight on the WNR shoulder ratio (0=single-peaked CAP, ->1=2nd peak)
maptol = gv('maptol', 167);  lc = gv('lc', [3.5 6.5]);
wshape = gv('wshape', 0);     % weight on abr_surface_obj's J = SHAPE-rms + b/c/d
                             % band hinges. DISTINCT from opts.surface, which uses the
                             % local surf_term (band violation in ms with a free Delta).
                             % Motivation (2026-07-25): the scalar d is NOISE-DOMINATED
                             % -- per-level slopes scatter 0.36-0.38, wider than the gap
                             % to target -- so fitting d fits noise. The SHAPE residual
                             % (rms log-resid vs tau=b*c^(-i)*f^(-d)) is ~0.132 even for
                             % the best model: large, well-determined, and unaddressed.
sobj   = gv('sobj', struct());% opts forwarded to abr_surface_obj (dband/cband/bband,
                             % wshape/wd/wc/wb) -- lets shape vs hinges be reweighted.
surface= gv('surface', 0);    % 1 = fit the RAW 16-point latency SURFACE (band loss, free Delta)
wsurf  = gv('wsurf', 1);      % weight on the surface band-violation RMS (ms)
dband  = gv('dband', [0 1.5]);% bounds (ms) on the free neural-delay offset Delta beyond the assumed 5 ms
wlcb   = gv('wlcb', 0.10);   % weight on the level_c band term (0 disables it)
tiptail= gv('tiptail', 0);   % 1 = fit the click-derived latency exponent d
dtarget= gv('dtarget', 0.40);% target d (1988/2013 band 0.39-0.41)
wd     = gv('wd', 1);        % weight on |d - dtarget|
hbm    = gv('hbmode', '');   % latency drive ('bm' = BM-peak). SET THIS for valid
                             % fits: the default drive rewards the WNR peak-switch
                             % artifact (that is what produced the bogus 0.413).
skipabr= gv('skipabr', 0);   % 1 = skip abr_metric in the START and END
                             % reconstructions. jointobj already skips it
                             % per-eval when no ABR term is active, but the two
                             % reconstructions called it unconditionally, and for
                             % m=4 it costs ~38 MIN EACH (16 tone-burst
                             % conditions). Set for a pure MAP fit. Leave 0 for
                             % tiptail fits: there the final abr_metric is the
                             % cross-check that exposed the linearized-d vs
                             % tone-burst-slope disconnect.
cheap  = gv('cheapstab', 0); % 1 = skip the PER-EVAL coupeig and rely on the
                             % tiptail CF-map screen as a gross-instability
                             % proxy. For m>=4, dof=3 gives N=8406 vs 2802 and
                             % dense eig scales ~N^3 (~27x), so per-eval coupeig
                             % is infeasible. Start/end coupeig still run.
                             % A finite click is NOT proof of stability, so this
                             % catches only GROSS failure -- the PINNED stability
                             % pair is what supplies the real static margin.
statol = gv('statol', 100);  % maxRe(all) tolerance. CALIBRATED 2026-07-25 from
                             % the champion baseline: HEALTHY models run to
                             % maxRe = +4.8 / +18.1 / +23.6 / +0.0, NOT the +0.0
                             % floor originally assumed -- a threshold of 1 would
                             % have rejected all four. Divergent models sit at
                             % +62747 / +102797 / +123306, so the separation is
                             % ~3 orders of magnitude and 100 is comfortably
                             % between. Still far tighter than the -40 margin
                             % used for the oscillatory (maxRe_osc) term.
pin    = gv('pin', struct());% fields forced onto the model and HELD OUT of the
                             % fit (e.g. m5o, chsz). Applied to base AND the warm
                             % start so pv0 carries them. Params absent from
                             % parnames (m5o/k5o/r5o) persist automatically;
                             % chsz(4) is index 34, so also drop it from fitidx.
out    = gv('out', sprintf('parfit26_nch%d.mat', nch));

base = modpar26(nch);
pa0  = base;
w = gv('warm', 'refit_c3_map.mat');
if (ischar(w)), if (exist(w,'file')), L=load(w); pa0=L.R.pa; end
elseif (isstruct(w)), pa0=w; end
% keep chamber count consistent with the chosen model
if (numel(pa0.chsz) ~= numel(base.chsz)), pa0=base; end
pf=fieldnames(pin);          % force pinned fields onto BOTH base and the warm
for i=1:numel(pf)            % start, so pv0 carries them and setpar_l cannot
    base.(pf{i})=pin.(pf{i});% overwrite them mid-fit
    pa0.(pf{i}) =pin.(pf{i});
end
if (~isempty(pf)), fprintf('  pinned: %s\n', strjoin(pf,', ')); end
pv0 = getpar_l(pa0);

chszderive = gv('chszderive', []); % index of the chsz element DERIVED from the
                             % others so sum(chsz)==2 identically. [] = off
                             % (chsz free, prior behaviour). Set to 3 to derive
                             % SV, the chamber SN's construction carves from.
                             % The derived index is removed from fitidx below --
                             % fitting it would be overwritten and is meaningless.
if (~isempty(chszderive))
    di = 30 + chszderive;    % getpar index of the derived chamber
    if (any(fitidx == di))
        fitidx = fitidx(fitidx ~= di);
        fprintf('  chsz(%d) DERIVED (sum=2 enforced); dropped index %d from fitidx\n', ...
                chszderive, di);
    else
        fprintf('  chsz(%d) DERIVED (sum=2 enforced)\n', chszderive);
    end
end
wstat  = gv('wstat', wcrit); % weight on the STATIC-divergence guard,
                             % wcrit*max(0, maxRe - statol). Defaults to wcrit
                             % so existing calls are unchanged. Separate because
                             % the two guards want different magnitudes: the osc
                             % term carries a large near-constant offset (it is
                             % live whenever maxRe_osc > -40) while the static
                             % term is zero until statol is exceeded, so one
                             % shared weight cannot size both. See the note at
                             % the J assembly.
P.chszderive=chszderive;
P.wgain=wgain; P.gainmin=gainmin; P.wslope=wslope; P.wmap=wmap; P.wcrit=wcrit; P.wstat=wstat; P.maptol=maptol; P.lclo=lc(1); P.lchi=lc(2); P.wlevel=wlevel;
P.wanchor=wanchor; P.anchor=anchor; P.wshoulder=wshoulder;
P.hbmode=hbm; P.tiptail=tiptail; P.dtarget=dtarget; P.wd=wd;
P.surface=surface; P.wsurf=wsurf; P.dlo=dband(1); P.dhi=dband(2); P.wlcb=wlcb;
P.wshape=wshape; P.sobj=sobj;
P.cheapstab=cheap; P.statol=statol;
obj = @(x) jointobj(x, pv0, fitidx, base, P);

abrstub = struct('slope',NaN,'level_c',NaN,'n_sub',NaN,'ok',true);
if (skipabr), m0 = abrstub; else, m0 = abr_metric(pa0, false); end
R0 = fdm26(struct('pa',pa0));
fprintf('parfit26 nch=%d  start: slope=%.3f  maperr=%.1f  level_c=%.2f  (fitting %d params)\n', ...
        nch, m0.slope, R0.maperr, m0.level_c, numel(fitidx));
if (skipabr), fprintf('  (skipabr=1: slope/level_c not computed -- pure MAP fit)\n'); end

ckpt = [out '.ckpt.mat'];   % checkpoint each iteration: two fits have died
                            % silently (machine sleep), losing all progress
op = optimset('Display','off','MaxFunEvals',maxfe,'MaxIter',maxfe, ...
              'OutputFcn',@(x,ov,st) outfun(x,ov,st,ckpt,pv0,fitidx));
t0=tic; [xopt,Jopt]=fminsearch(obj, pv0(fitidx), op); wall=toc(t0);

pv=pv0; pv(fitidx)=xopt; pa=setpar_l(base,pv,chszderive);
if (~isempty(hbm)), pa.hbmode=hbm; end
if (skipabr), mf=abrstub; else, mf=abr_metric(pa,false); end
evalc('S=tdm26(''coupeig'',struct(''pa'',pa));'); Rf=fdm26(struct('pa',pa));
fprintf('\n=== parfit26 nch=%d result (%.0f s, J final %.4f) ===\n', nch, wall, Jopt);
fprintf('  ABR slope d : %.3f  (target 0.413)\n', mf.slope);
fprintf('  maperr      : %.1f  (target <= %.0f)\n', Rf.maperr, maptol);
fprintf('  level_c     : %.2f  (target 5)\n', mf.level_c);
% STABILITY LINE. Reports the FULL maxRe and bases the verdict on it, with
% maxRe_osc shown alongside as secondary. It used to print maxRe_osc ALONE and
% label it sub-critical/UNSTABLE from that value, which is actively misleading:
% maxRe_osc filters to oscillatory modes above 1 kHz and is BLIND to static
% (f=0) divergence. On refit_m4_map (2026-07-28) it printed "osc-maxRe -2.3
% (sub-critical)" for a model whose true maxRe was +27.6 -- unstable, and the
% result void. The same blindness once reported -486 for a model whose time
% march reached 3e296. maxRe(all) is the arbiter; osc is a diagnostic.
% Threshold 24 matches score26's healthy band so the two agree.
vst = {'UNSTABLE', 'healthy'};
fprintf('  maxRe       : %+.1f (%s)  [osc %+.1f], n_sub=%d%s\n', ...
        S.maxRe, vst{1 + (S.maxRe < 24)}, S.maxRe_osc, mf.n_sub, ...
        tern(S.maxRe >= 24 && S.maxRe_osc < 0, '   <== osc alone would say sub-critical', ''));
if (surface)
    [st,dbest,nv]=surf_term(mf,P);
    fprintf('  SURFACE     : band-RMS %.3f ms over %d/%d cells\n', st, nv, numel(mf.lat));
    fprintf('  Delta (pred): %.2f ms  =>  neural delay a = %.2f ms  (2013 assumed 5.00)\n', dbest, 5+dbest);
    R.delta=dbest; R.surf_rms=st; R.surf_nvalid=nv;
end
R.pa=pa; R.pv=pv; R.mf=mf; R.S=S; R.Rf=Rf; R.fitidx=fitidx; R.nch=nch; R.opts=opts;
save(out,'R'); fprintf('saved %s\n', out);
end

% ---- joint objective: tuning (fdm26) + slope (abr_metric) + stability (coupeig) ----
function J = jointobj(x, pv0, fitidx, base, P)
pv=pv0; pv(fitidx)=x; pa=setpar_l(base,pv,P.chszderive);
if (~isempty(P.hbmode)), pa.hbmode=P.hbmode; end   % setpar_l rebuilds pa from
                                                   % modpar26, dropping extras
try
    % The time-domain ABR metric costs ~150 s per evaluation, so skip it entirely
    % when NO ABR-dependent term is active (e.g. a pure tuning/MAP fit, Test A).
    % fdm26 + coupeig alone are seconds, making such a fit ~100x cheaper.
    tt=[];
    if (P.tiptail)
        % Click-derived exponent d: ~11 s/eval vs ~100 s for the tone-burst
        % metric, so this supports a far larger budget. Skips abr_metric.
        tt=tiptail_metric(pa);
        if (~tt.ok), J=1e6; return; end
        needm=false;
    else
        needm = P.surface || P.wshape~=0 || P.wslope~=0 || P.wlevel~=0 || P.wanchor~=0 || ...
                P.wshoulder~=0 || P.wlcb~=0;
    end
    m=[];
    if (needm)
        m=abr_metric(pa,false);              % TDM ABR slope + level_c (reliable)
        if (~m.ok), J=1e6; return; end
    end
    if (P.cheapstab)
        % Per-eval coupeig skipped (see opts.cheapstab). tiptail_metric's CF-map
        % screen + R2 gate already rejects grossly unstable/degenerate models via
        % the tt.ok test above, which is the proxy relied on here. Neutral values
        % so the stability terms below contribute nothing.
        S.maxRe_osc = -Inf; S.maxRe = 0;
    else
        evalc('S=tdm26(''coupeig'',struct(''pa'',pa));');   % sub-criticality
    end
    % CONSOLIDATED (2026-07-29). score26('fast') ALREADY runs fdm26 and coupeig
    % internally, so calling it for amp_gain on top of the separate calls above
    % paid for each TWICE. Measured consequence: a 20-evaluation m=4 fit had not
    % cleared setup after 2h49m. One call now supplies maperr, maxRe, maxRe_osc
    % AND amp_gain, so the amplifier guard costs only the click block (two extra
    % tdm26 runs) rather than doubling the whole evaluation.
    agn = NaN;
    if (P.wgain ~= 0)
        Sg = [];
        try, Sg = score26(pa,'fast',false); catch, end
        if (isstruct(Sg) && isfield(Sg,'maperr'))
            Rf = struct('maperr', Sg.maperr);
            if (isfield(Sg,'maxRe') && isfinite(Sg.maxRe)), S.maxRe = Sg.maxRe; end
            if (isfield(Sg,'maxRe_osc') && isfinite(Sg.maxRe_osc)), S.maxRe_osc = Sg.maxRe_osc; end
            if (isfield(Sg,'amp_gain') && isnumeric(Sg.amp_gain) && ~isempty(Sg.amp_gain))
                agn = double(Sg.amp_gain(1));
            end
        else
            J = 1e6; return;                 % score26 failed: reject this point
        end
    else
        Rf=fdm26(struct('pa',pa));           % fd tuning/MAP error (delegated)
    end
    Jm=0;                                    % all ABR-dependent terms
    if (P.tiptail)
        Jm = P.wd*abs(tt.d - P.dtarget) ...
           + 0.02*max(0, 12-tt.chi) ...      % keep a distinct tip at high CF
           + 0.10*max(0, 4-tt.nvalid);       % keep the CF map intact. Threshold 4,
                                             % not 5: tiptail_metric already fails
                                             % below 4, and the 0.6-16 kHz window
                                             % legitimately admits only 4 places
                                             % for the native 2-chamber -- a 5
                                             % threshold penalized that filter
                                             % artifact more than |d-target|.
    elseif (needm)
    fi=find(abs(m.f-P.anchor(1))<0.01,1); li=find(abs(m.slv-P.anchor(2))<0.1,1);  % anchor grid point
    la=m.lat(fi,li); if (~isfinite(la)), la=m.abr(fi,li)+5; end   % sub-threshold at anchor -> 5 ms penalty
    sh=m.shoulder; if (~isfinite(sh)), sh=1; end                  % WNR shoulder ratio (0=single-peaked)
    if (P.surface)
        % REFORMULATED OBJECTIVE. Fit the RAW 16-point latency surface to
        % tau_data(f,i)+Delta, with Delta -- the neural delay IN EXCESS of the
        % assumed 5 ms -- FREE and bounded, under a band loss spanning the 1988
        % and 2013 power laws. Rationale: the 2013 paper's 5 ms (1 ms synaptic +
        % 4 ms wave I->V) is only established as LEVEL-independent; its
        % FREQUENCY-independence is an explicit simplifying assumption. Forcing
        % the model to absorb it biases the slope. Freeing Delta makes the neural
        % delay a model PREDICTION, and the surface uses all 16 cells rather than
        % compressing them into two summary numbers (slope, level_c).
        Jm = P.wsurf    *surf_term(m,P) ...
           + P.wshoulder*sh;
    else
        Jm = P.wslope *abs(m.slope-0.413) ...
           + P.wlevel *abs(100*(m.level_c^(1/100)-1) - 1.62) ... % level dependence as %/dB (target 1.62)
           + P.wanchor*abs(la - m.abr(fi,li)) ...                % absolute-latency anchor at (freq,level), ms
           + P.wshoulder*sh ...                                  % single-peaked WNR (no 2nd-DOF shoulder)
           + P.wlcb   *(max(0, P.lclo - m.level_c) + max(0, m.level_c - P.lchi));
    end
    % SHAPE term: abr_surface_obj fits all 16 cells JOINTLY to tau=b*c^(-i)*f^(-d)
    % and returns J = shape-rms + hinge losses on b,c,d that are ZERO inside the
    % published uncertainty. Added on top of whichever branch ran, so it can be
    % used alone (all other weights 0) or as an extra constraint.
    if (P.wshape ~= 0)
        [Js,~] = abr_surface_obj(m, P.sobj);
        if (~isfinite(Js)), J=1e6; return; end
        Jm = Jm + P.wshape * Js;
    end
    end
    % AMPLIFIER SHORTFALL -- see the note at the top. Skipped entirely when
    % wgain==0, so callers who opt out do not pay for the extra score26 call.
    Jg = 0;
    if (P.wgain ~= 0)
        if (~isfinite(agn)), J = 1e6; return; end  % no measurable amplifier
        Jg = P.wgain * max(0, P.gainmin - agn);
    end
    J = Jm ...
      + P.wcrit *max(0, S.maxRe_osc + 40) ...
      + P.wstat *max(0, S.maxRe - P.statol) ...
      + P.wmap  *max(0, Rf.maperr - P.maptol) ...
      + Jg;
    % SEPARATE WEIGHTS (2026-07-28). These two guards previously SHARED wcrit,
    % and that sharing forced a bad trade. The osc term is nonzero whenever
    % maxRe_osc > -40, so at a typical osc of -2.5 it contributes wcrit*37.5 --
    % a large, nearly constant offset. The static term only bites above statol.
    % Raising wcrit enough for the static guard to matter therefore inflates the
    % osc term until it dominates the objective.
    % Measured in refit_m4_map2 (wcrit=0.05, statol=20): the final J of 2.2123
    % decomposed as map 0.337, osc 1.875, static 0.000 -- the osc guard was 5.6x
    % the metric being fitted, so most of the search went into pushing maxRe_osc
    % further sub-critical rather than improving the map.
    % With wstat separate, the static guard can be made to bite (it needs to
    % outweigh the map gain it would otherwise buy: a 46-point maperr gain is
    % worth wmap*46, so wstat*(maxRe-statol) must exceed that at the maxRe being
    % prevented) WITHOUT inflating the osc term.
    % wstat defaults to wcrit, so any existing call reproduces its old behaviour.
    % The second term is the STATIC-divergence guard. maxRe_osc filters to
    % oscillatory modes above 1 kHz and is BLIND to f=0 divergence: m=4 at
    % gam=0.5 reported maxRe_osc=-486 ("stable") while maxRe=+16862 and the
    % time march blew up to 3e296. Stable configs sit at maxRe=+0.0 (spurious
    % perturbed-zero boundary modes), hence the small positive P.statol.
    if (~isfinite(J)), J=1e6; end
catch
    J=1e6;
end
end

function stop=outfun(x,ov,state,ckpt,pv0,fitidx)
stop=false;
if (strcmp(state,'iter'))
    fprintf('  iter %2d: J=%.4f  (fevals %d)\n', ov.iteration, ov.fval, ov.funccount);
    try   % best-effort checkpoint; never let a save failure stop the fit
        C.x=x; C.J=ov.fval; C.iter=ov.iteration; C.fevals=ov.funccount;
        C.pv0=pv0; C.fitidx=fitidx; C.when=datestr(now); %#ok<TNOW1,DATST>
        save(ckpt,'C');
    catch
    end
end
end

% ---- target band for the forward-latency surface --------------------------
% Published power laws tau = b*c^(-i)*f^(-d), i = dB/100, f = kHz:
%   Neely et al. 1988      b=12.90 c=5.00 d=0.413
%   Rasetshwane et al 2013 b=12.63 c=5.34 d=0.390   (replication)
% The two bracket the experimental uncertainty, so we penalize only excursions
% OUTSIDE the band they span -- a hinge/band loss rather than a point target.
function [tlo,thi]=abr_target_band(f,slv)
[F,I]=ndgrid(f(:), slv(:)/100);
t88 = 12.90*(5.00.^(-I)).*(F.^(-0.413));
t13 = 12.63*(5.34.^(-I)).*(F.^(-0.390));
tlo = min(t88,t13); thi = max(t88,t13);
end

% ---- surface band-violation with a FREE neural-delay offset Delta ---------
% Returns the RMS band violation (ms) of (model latency + Delta), minimized over
% Delta in [P.dlo P.dhi]. Delta is the neural delay IN EXCESS of the assumed
% 5 ms; the published forward latencies are wave-V minus a CONSTANT 5 ms whose
% frequency-independence is an assumption, so Delta absorbs exactly that bias
% and is reported as a model prediction (a = 5 + Delta).
function [term,dbest,nv]=surf_term(m,P)
[tlo,thi]=abr_target_band(m.f, m.slv);
L=m.lat; ok=isfinite(L) & L>0;          % drop sub-threshold / failed-detector cells
nv=sum(ok(:));
if (nv < 8), term=1e3; dbest=NaN; return; end
ds=linspace(P.dlo,P.dhi,301); term=inf; dbest=P.dlo;
for k=1:numel(ds)
    T=L+ds(k);
    v=max(0, tlo-T) + max(0, T-thi);    % 0 inside the band, distance to nearest edge outside
    r=sqrt(mean(v(ok).^2));
    if (r<term), term=r; dbest=ds(k); end
end
end

% ---- 30 impedance params (getpar order) + chamber sizes (the 3-chamber lever) ----
function nm=parnames()
nm={'k1o','r1o','m1o','k2o','r2o','m2o','k3o','r3o','k4o','aco', ...
    'k1e','r1e','m1e','k2e','r2e','m2e','k3e','r3e','k4e','ace', ...
    'k1q','r1q','m1q','k2q','r2q','m2q','k3q','r3q','k4q','acq'};
end
function pv=getpar_l(pa)
nm=parnames(); nc=numel(pa.chsz); pv=zeros(1,30+nc);
for i=1:30, pv(i)=pa.(nm{i}); end
pv(31:30+nc)=pa.chsz(:)';
end
function pa=setpar_l(pa,pv,dvi)
% dvi = index of the chsz element DERIVED from the others so that
% sum(chsz) == 2 identically. Empty or absent leaves chsz free (prior
% behaviour). This enforces the constraint INSIDE the parameter mapping, so no
% search step can violate it -- a post-hoc check or a renormalisation could not
% do that, and renormalisation was explicitly rejected earlier because it makes
% added chambers resize existing ones.
% SN's construction: chambers are carved from SV, so SV is the natural derived
% element (index 3). Both readings of the dual bisection are reachable --
% "both cuts in SV" gives [1, s, 1-s-c, c] and "one cut in ST, one in SV" gives
% [1-s, s, 1-c, c] -- and with only sum=2 enforced the fit can land anywhere
% between them, which lets the result choose the construction rather than
% assuming it.
nm=parnames(); nc=numel(pa.chsz);
for i=1:30, pa.(nm{i})=pv(i); end
pa.chsz=pv(31:30+nc);
if (nargin>=3 && ~isempty(dvi) && dvi>=1 && dvi<=nc)
    o = true(1,nc); o(dvi) = false;
    pa.chsz(dvi) = 2 - sum(pa.chsz(o));
end
end
function v=subsref_default(s,f,d), if (isfield(s,f)), v=s.(f); else, v=d; end, end
function s=tern(c,a,b), if c, s=a; else, s=b; end, end
