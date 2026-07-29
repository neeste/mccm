function S = tdm26(stcfg,nch,dsp1,dsp2) % gem
% Multi-chamber cochlear-model responses
%
%   USAGE: S = tdm26(STCFG,NCH,DSP1,DSP2)
%
%   NCH specifies the number of fluid chambers
%       When NCH=0, computes only the middle-ear stage.
%   STCCFG specifies the type of stimulus, for example:
%       0        - selects a click,
%       F        - selects a tone with frequency = F kHz,
%     'file.wav' - selects a waveform stimulus.
%   Otherwise, STCFG may select an experimental protocol:
%     'tbabr'    - whole-nerve-response (WNR) latency to tonebursts,
%     'fwdmsk'   - tone-on-tone forward masking of WNR.
%     'ldngrw'   - loundess growth of single tone,
%     'dpoae'    - distortion-product otoacoustic emission (NYI).
%     'param'    - replace model parameters
%     'ecochg'   - ECochG simulation
%   DSP1 = upper subplot [1=Vr 2=Pe 3=Ps 4=WNR]
%   DSP2 = lower subplot [1=HB 2=NR 3=LD 4=HC]

if (nargin<1), stcfg = 0; end
if (nargin<2), nch = 0;   end
if (nargin<3), dsp1 = 0;  end
if (nargin<4), dsp2 = 0;  end

% process protocol
ts = check_protocol(stcfg,nch,dsp1,dsp2);
if (ts.done), if (nargout), S = ts; end; return; end

% set default display
if (nargin<3), dsp1 = 1;  end
if (nargin<4), dsp2 = 2;  end

% Allow a parameter STRUCT in the NCH slot for non-protocol stimuli (click,
% tone), mirroring what the protocols already accept. This lets a response be
% run with non-default parameters -- e.g. m=2 driven by the cel26c3 impedance
% set -- which separates CHAMBER COUNT from PARAMETER SET.
pa_in=[]; if (isstruct(nch)), pa_in=nch; nch=pa_in.m; end

% compute time-domain model response
if (isempty(pa_in))
    [pa,sav,cur,ts] = tdm_init(stcfg,nch,dsp1,dsp2,ts);
else
    [pa,sav,cur,ts] = tdm_init(stcfg,nch,dsp1,dsp2,ts,pa_in);
end
[sav,cur,ts,cp] = tdm_step(pa,sav,cur,ts);
       [mbf,bf] = tdm_show(pa,sav,cur,ts,cp);

if (nargout)
    [~,~,pe,ps,~,~,ve]=fetch_sav(1,sav);
    sav.gme=ps./pe; sav.zme=pe./ve;
    sav.pa = pa; sav.ts = ts; sav.cp = cp;
    sav.bf = bf; sav.mbf = mbf; S = sav;
end
end

%==========================================================
% CORE SIMULATION FUNCTIONS
%==========================================================

function [pa,sav,cur,ts]=tdm_init(stcfg,nch,dsp1,dsp2,ts,pa)
if (nargin<6), pa=modpar26(nch); end
if (nch<1), pa.ihcv = 0; end
if (pa.nmev==1 && pa.stim==3), pa.stim=2; fprintf('moved earphone stimulus to eardrum\n'); end

dof = resolve_dof(pa);   % Stage 3: DOF count is a PARAMETER (see resolve_dof)
% m=3b (pa.d3int) has THREE DOFs on THREE chambers: d3 is INTERNAL, a local
% mass-spring with no fluid compartment, exactly as d2 is in m=2. It is the
% missing rung between m=3 and m=4.
pa.dof = dof; pa.ncp = pa.n * dof;
if (pa.nmev<1), pa.nmev = 1; end
nsv = pa.ncp + pa.nmev;
cur.d = zeros(nsv,1); cur.v = zeros(nsv,1);
cur.hb = zeros(pa.n,1); cur.hc = zeros(pa.n,1);
cur.nr = zeros(pa.n,1); cur.ld = zeros(pa.ldne,1);
cur.stp = 0; cur.cyc = 0; cur.wnr = 0; cur.ohcp = 0; cur.ohcbm = 0;
cur.vst = 0; cur.pst = 0; cur.ped = 0;
cur.ved = 0; cur.vep = 0; cur.stm = 0;
cur.nrme = zeros(pa.nmev,1); cur.qme = zeros(4,1);

nt = pa.ntsf;
sav.vep = zeros(nt,1); sav.ped = zeros(nt,1); sav.ved = zeros(nt,1);
sav.pst = zeros(nt,1); sav.vst = zeros(nt,1); sav.wnr = zeros(nt,1);
sav.ohcp = zeros(nt,1);
sav.ohcbm = zeros(nt,1);   % OHC work against BM velocity
sav.d1 = zeros(nt,length(pa.isv)); sav.d2 = zeros(nt,length(pa.isv));
sav.d3 = zeros(nt,length(pa.isv));   % 3rd DOF (OC height) -- 4-chamber
sav.hb = zeros(nt,pa.n); sav.hc = zeros(nt,pa.n);
sav.nr = zeros(nt,pa.n); sav.ld = zeros(nt,pa.ldne);
sav.nch = nch; sav.isv = pa.isv;

ts=ts_init(pa,sav,cur,ts,1,2,stcfg,dsp1,dsp2);
if (ts.tt == 0), return; end
ts.smx = 0;
nw = ts.nw; nt = ts.nt / nw; dt = ts.dt * nw / 1000;
sav.f = (0:((nt/2) - 1))' / (1000 * dt * nt);
sav.smx = ts.smx;
end

function ts=ts_init(pa,sav,cur,ts,tsfig,svfig,stcfg,dsp1,dsp2)
if (~dsp1 && ~dsp2), tsfig = 0; end
ts.tsdsp = tsfig;
ts.svdsp = svfig;
ts.stcfg = stcfg;
ts.lab = sprintf('nch=%d',pa.m);
[pa,stwav,tt,nt,smx]=stim_init(pa,stcfg);
ts.stwav = stwav;
ts.dsp1 = dsp1;
ts.dsp2 = dsp2;
ts = plot_ts(pa,sav,cur,ts,tt,smx);
ts.tt = tt;
ts.nt = nt;
ts.dt = tt / nt;
ts.nw = pa.ntsw;
ts.ntsf = floor(ts.nt/ts.nw);
end

function [pa,stwav,tt,nt,smx]=stim_init(pa,stcfg)
sr=1/pa.dt; dt=2e-5; nt=pa.nstp; nn=round(nt*pa.dt/dt); tt=nn*dt*1000;
if (iscell(stcfg)), [yy,nt,tt]=waveform(stcfg,sr);
elseif (ischar(stcfg)), [yy,nt,tt]=waveform(stcfg,sr);
elseif (stcfg==0), n1=1+round(0.001/dt); n2=1+round(0.002/dt); yy=short_chirp(nn,n1,n2,nt);
elseif (size(stcfg,1)<=2), [yy,nt,tt]=multi_tone(stcfg,tt,nt,0);
elseif (stcfg(1,1)>0), [yy,nt,tt]=multi_toneburst(stcfg,tt,nt);
else, st=stcfg(4,1); rc=stcfg(5,1); [yy,nt,tt]=band_limit(stcfg,tt,nt,st,rc);
end
stwav=yy * pa.stgain; smx=max(max(abs(stwav)),1e-9);
end

function [sav,cur,ts,cp]=tdm_step(pa,sav,cur,ts)
cp=[];
if (pa.m<0), ts=plot_spec(pa,sav,ts); end
if (ts.done), return; end

me=midear(pa);
[cp,me]=macro26(pa,me);
[cur,nf]=neural(pa,cur);

nxt = applst(pa,ts,cur);
sav = savest(pa,sav,cur,ts);

while (cur.cyc < pa.ncyc)
    nxt = inctim(cur,nxt,ts);
  	[cur,ac] = accel(pa,cp,me,cur,cur);
    cur = hc_step(pa,cur,nf);
    sav = savest(pa,sav,cur,ts);
  	nxt = update(pa,cur,nxt,ac,ac);
    nxt = applst(pa,ts,nxt);
  	for imp=1:pa.nimp
		[cur,an] = accel(pa,cp,me,cur,nxt);
     	nxt = update(pa,cur,nxt,ac,an);
    end
	cur=nxt;
end
end

function st=applst(pa,ts,st)
if (isempty(ts.stwav)), return; end
if (st.stp>=ts.nt), stm = 0; else, stm = ts.stwav(st.stp+1); end
switch (pa.stim)
    case 1, st.pst=stm;
    case 2, st.ped=stm;
    case 3, st.vep=stm;
end
st.stm=stm;
st.qme(pa.stim)=stm;
end

function nxt=inctim(cur,nxt,ts)
nxt.stp = cur.stp + 1;
nxt.cyc = cur.cyc;
if (nxt.stp >= ts.nt )
   nxt.stp = nxt.stp - ts.nt;
   nxt.cyc = nxt.cyc + 1;
end
end

function [cur,a]=accel(pa,cp,me,cur,st)
a = zeros(size(st.d)); stm = st.stm;
nme = pa.nmev; ime = pa.ncp + (1:nme);
d0 = st.d(ime); v0 = st.v(ime);
ist = me.idx(2);
[ss,ii,gam,ohcp,ohcbm]=micro26(pa,cp,st);

s0 = stm * me.smes - (me.smek * d0 + me.smer * v0);
prw = (pa.rrw * v0(ist) + pa.krw * d0(ist)) * pa.arw;
m = pa.m; n = pa.n;

% Bi-directional (reciprocal) stapes port -- NOW THE DEFAULT (pa.me_recip=0
% restores the legacy asymmetric port). The SAME coefficient sqrt(alfx) is used
% BOTH ways, so the port is its own transpose while the round-trip product
% afw*arv = alfx is PRESERVED (loop gain, hence stability, unchanged) and the
% emission couples outward reciprocally with the drive.
% LEGACY (me_recip=0) drove the cochlea with cp.alfx but loaded the stapes with
% 1 -- the only non-reciprocal coupling in the model; every other port uses the
% symmetric D/D^T pair xpnd_q/fold_p. NB: symmetrizing by scaling ONLY the
% reverse path to alfx raises the loop gain to alfx^2 and DIVERGES; the even
% sqrt split is what keeps the round trip (and stability) invariant.
recip=1; if (isfield(pa,'me_recip')), recip=pa.me_recip; end   % DEFAULT: reciprocal
if (recip), afw=sqrt(cp.alfx); arv=afw; else, afw=cp.alfx; arv=1; end

if (pa.stim==1), qst = stm;
elseif (m), qst = -afw * (s0(ist) + prw);
else, qst = me.qstr * v0 + me.qstk * d0;
end

if (m), qq=xpnd_q(qst,ss,m,n,cp); pp = cp.aa \ qq;
else, pp = qst; end

% CALCULATE DIFFERENTIAL PRESSURE FOR THE STAPES
if (m>=3)
    p_dif = pp(1) - pp(3); % P_ST - P_SV (chamber order 1=ST,2=SS,3=SV[,4=CL])
else
    p_dif = pp(1);         % 1 and 2 chamber already use P_diff
end

% Bi-directional (reciprocal) stapes coupling: the forward path injects the ME
% state into the cochlea with coefficient cp.alfx (qst = -alfx*(s0+prw), above),
% so the reverse path -- the base differential pressure loading the stapes --
% must be its TRANSPOSE and carry the SAME coefficient. With rcp=1 (default,
% legacy) the stapes port is the only non-reciprocal coupling in the model
% (every other port uses the symmetric D/D^T pair xpnd_q/fold_p), which
% under-couples the emission outward and mis-terminates the base.
s0(ist) =  s0(ist) + arv*p_dif + prw;

if me.use_decomp, a0 = me.mass_decomp \ s0; else, a0 = me.smem \ s0; end

a(ime)=a0;
if (m), a = fold_p(pp,m,n,a,ss,ii,cp);
else, a(ii(:,1)) = -a0(ist); end

if (pa.stim==2), cur.vep = cur.ped;
else, cur.ped = me.pedm * a0 + me.pedr * v0 + me.pedk * d0; end

cur.ohcp = ohcp; cur.ohcbm = ohcbm;
cur.ved = me.ved * v0; cur.vst = me.vst * v0;
cur.pst = pp(1); cur.gam = gam; st.a = a;
end

function nxt=update(pa,cur,nxt,cur_a,nxt_a)
dt2 = pa.dt / 2;
nxt.v = cur.v + (cur_a + nxt_a) * dt2;
nxt.d = cur.d + (cur.v + nxt.v) * dt2;
nxt.vi = cur.vi; nxt.ss = cur.ss; nxt.nr = cur.nr; nxt.ld = cur.ld;
end

function sav=savest(pa,sav,cur,ts)
istp = cur.stp; n=pa.n; nw=ts.nw; kw=mod(istp,nw);
if (kw == 0)
   i = 1 + mod(round(istp / nw), ts.nt);
   sav.vep(i) = cur.vep; sav.ped(i) = cur.ped; sav.ved(i) = cur.ved;
   sav.ohcp(i) = cur.ohcp; sav.ohcbm(i) = cur.ohcbm;
   sav.pst(i) = cur.pst; sav.vst(i) = cur.vst; sav.wnr(i) = cur.wnr;
   sav.d1(i,:) = cur.d(pa.isv).'; sav.d2(i,:) = cur.d(pa.isv+n).';
   if (pa.dof>=3), sav.d3(i,:) = cur.d(pa.isv+2*n).'; end   % m>=4 OR m=3b:
                        % gate on the DOF COUNT, not the chamber count. Gating
                        % on pa.m>=4 left m=3b computing d3 correctly but never
                        % SAVING it, so it read as all zeros -- a save-guard
                        % artifact that looks exactly like 'd3 is not driven'.
   sav.hb(i,:) = cur.hb; sav.hc(i,:) = cur.hc;
   sav.nr(i,:) = cur.nr; sav.ld(i,:) = cur.ld;
   ts.smx = sav.smx;
   ts = plot_ts(pa,sav,cur,ts,i);
   sav.smx = ts.smx;
end
end

%==========================================================
% PHYSICS MODULES
%==========================================================

function me=midear(pa)
nme=pa.nmev;
me.pedm=zeros(1,nme); me.pedr=zeros(1,nme); me.pedk=zeros(1,nme);
me.smes=zeros(nme,1); me.qstr=zeros(1,nme); me.qstk=zeros(1,nme);
me.smem=zeros(nme,nme); me.smer=zeros(nme,nme); me.smek=zeros(nme,nme);

if (pa.nmev==1)
    me.smem = 1;
    me.smer = (pa.rme / pa.mme);
    me.smek = (pa.kme / pa.mme);
    me.qstr = pa.rme; me.qstk = pa.kme; me.smes = -1;
    me.ied = 1; me.ist=1; me.idx = [1 1]; me.vst = [0 0 0 1];
    return
end

iep=1; icp=2; ima=3; ist=4;
me.idx=[icp ist];
cr = pa.cep; ad = pa.adi; ac = pa.acp; au = pa.aed; am = pa.ama; as = pa.ast;
gm = pa.gme; md = pa.mdi; mc = pa.mcp; mu = pa.med; mm = pa.mma; ms = pa.mst;
ru = pa.red; rm = pa.rma; ri = pa.rim; rs = pa.rst; rc = pa.rcp; ra = pa.rfz;
rd = pa.rdi + pa.rep; ku = pa.ked; km = pa.kma; ki = pa.kim; ks = pa.kst;
ka = pa.kcp; kd = pa.kdi;

me.smem(iep,iep) = md / ad^2;
me.smem(icp,icp) = mc / ac^2 + mu  / au^2;
me.smem(icp,ima) = -mu  / au^2; me.smem(ima,icp) = -mu  / au^2;
me.smem(ima,ima) = mm / am^2 + mu  / au^2;
me.smem(ist,ist) = ms / as^2;

me.smer(iep,iep) = (rd + ra) / ad^2;
me.smer(iep,icp) = -ra / (ad * ac); me.smer(icp,iep) = -ra / (ad * ac);
me.smer(icp,icp) = (rc + ra) / ac^2 + ru / au^2;
me.smer(icp,ima) = -ru / au^2; me.smer(icp,ima) = -ru / au^2;
me.smer(ima,icp) = -ru / au^2;
me.smer(ima,ima) = (rm +  gm^2 * ri) / am^2 + ru / au^2;
me.smer(ima,ist) = -gm * ri / (am * as); me.smer(ist,ima) = -gm * ri / (am * as);
me.smer(ist,ist) = (rs + ri) / as^2;

me.smek(iep,iep) = (kd + ka) / ad^2;
me.smek(iep,icp) = -ka / (ad * ac); me.smek(icp,iep) = -ka / (ad * ac);
me.smek(icp,icp) = ka / ac^2 + ku / au^2;
me.smek(icp,ima) = -ku / au^2; me.smek(ima,icp) = -ku / au^2;
me.smek(ima,icp) = -ku / au^2;
me.smek(ima,ima) = (km +  gm^2 * ki) / am^2 + ku / au^2;
me.smek(ima,ist) = -(gm * ki) / (am * as); me.smek(ist,ima) = -(gm * ki) / (am * as);
me.smek(ist,ist) = (ks + ki) / as^2;

try me.mass_decomp = decomposition(me.smem); me.use_decomp = true;
catch, me.use_decomp = false; end

if (pa.stim==3), me.smes(iep) = -(cr / ad);
elseif (pa.stim==2), me.smes(icp) = -1; me.smes(ima) = -1; end
me.pedm(icp) = mu / au^2; me.pedr(icp) = ru / au^2; me.pedk(icp) = ku / au^2;
me.pedm(ima) = -mu / au^2; me.pedr(ima) = -ru / au^2; me.pedk(ima) = -ku / au^2;
qstrst = (rs + ri) / ms; qstkst = (ks + ki) / ms;
qstrma = (gm * as * ri) / (am * ms); qstkma = (gm * as * ki) / (am * ms);
if (pa.m)
    mst = ms / as^2;
    me.qstr(ist) = -mst * qstrst; me.qstk(ist) = -mst * qstkst;
    me.qstr(ima) =  mst * qstrma; me.qstk(ima) =  mst * qstkma;
else
    asco = as + pa.mco * as^2 / ms;
    me.qstr(ist) = (pa.rco - pa.mco * qstrst) / asco;
    me.qstk(ist) = (       - pa.mco * qstkst) / asco;
    me.qstr(ima) = (         pa.mco * qstrma) / asco;
    me.qstk(ima) = (         pa.mco * qstkma) / asco;
    me.mem = 0;
end
me.mst = (ms / as);
me.ved = [0 pa.ME_VEP_SCALE 0 0];
if (pa.m==1), me.vst = [0 0 0 pa.ME_VST_SCALE_PHANTOM]; else, me.vst = [0 0 0 1]; end
end

% cochlea/imped/force_cp EXTRACTED to macro26.m / imped26.m / micro26.m
% (Stage 4). tdm26 keeps the time march and the solver plumbing.

function [cur,nf]=neural(pa,cur)
n = pa.n;
cur.vohc = zeros(n,1); cur.dvohc = zeros(n,1); cur.vi = zeros(n,1);
cur.cc = zeros(n,1); cur.ss  = ones(n,1); cur.nr  = zeros(n,1);
nf.nihc = pa.IHC_SUBSTEPS; nf.sf = pa.ihcsf;
dt = pa.dt * nf.nihc;
nf.kk = exp(-dt / pa.ihctc);
nf.rr = pa.ihcrr * (dt / 2e-6)^(1/pa.ihcex);
nf.dr = pa.ihcdr * (dt / 2e-6)^(1/pa.ihcex);
ehw=floor(pa.n/(pa.xl*2/pa.ew)); i1k = pa.isv(3);
nf.ens = (i1k - ehw) : (i1k + ehw);
end

function cur=hc_step(pa,cur,nf)
n = pa.n; ii=1:n;
dd = cur.d(ii) - cur.d(ii+n);
met=pa.met; hco=pa.hco; hcs=pa.hcs;
b0=1./(1+exp(hco/hcs)); b1=1./(1+exp((hco-dd)/hcs));
cur.hc=(b1-b0)*hcs*met;
v1 = cur.v(ii); v2 = cur.v(ii+n);        % DOF-1 (BM) and DOF-2 (2nd radial) velocities
if (pa.ihcv==0), cur.hb = dd;
elseif (isfield(pa,'hbmode') && strcmp(pa.hbmode,'bm')),   cur.hb = v1;    % diag: DOF-1 (BM) drive only
elseif (isfield(pa,'hbmode') && strcmp(pa.hbmode,'dof2')), cur.hb = -v2;   % diag: DOF-2 drive only
else, a=1; if (isfield(pa,'hbalpha')), a=pa.hbalpha; end; cur.hb = v1-a*v2; % shear w/ 2nd-DOF weight a (dflt 1)
end
cur.nr = cur.hb.*pa.synpro;
if (pa.ihceq<1), return; end
if ((nf.nihc>1) && mod(cur.stp,nf.nihc)), return; end
dt = pa.dt * nf.nihc;
s = cur.hb; s = max(0,s); s = (1 - nf.kk) .* s + nf.kk .* cur.vi; cur.vi = s;
vi = s * nf.sf; ss = cur.ss; nr = cur.nr;
ss = ss + (nf.rr * (1 - ss) - nf.dr * nr) * dt; ss = max(0, min(ss, 1));
nr = ss.^pa.ihcex .* vi;
cur.ss = ss; cur.nr = nr; cur.wnr = sum(cur.nr) * pa.nrgn;
ew = pa.ldew/(pa.xl*10); ne = pa.ldne; sc = pa.ldsc;
cur.ld = soft_max(loudness(cur.nr,ew,ne,sc).^(1/3),10);
end

function x=soft_max(x,m)
x=x./(m+sum(x));
end

function [pld,enr]=loudness(nr,ewn,ne,c)
nx=size(nr,1); kw=round(nx*ewn)-1; enr=zeros(ne,1); p=1/3;
for k=1:ne, k1=1+round((k-1)*(nx-1)/(ne-1)); k2=min(k1+kw,nx); kk=k1:k2; enr(k)=mean(mean(abs(nr(kk,:)).^p)); end
N=enr.^(1/p); pld=(N.^c(1)./(1+N.^c(2)/c(3)))*c(4);
end

function qx=xpnd_q(qst,ss,m,n,cp)
% MACRO/MICRO INTERFACE (Stage 2). This was ~45 lines of hand-written
% per-chamber volume injections branching on m; it is one matrix operation with
% the topology supplied by macro_couple:
%
%     q(:,c) = sum_k Dq(c,k) * mu(:,k) .* s(:,k)  +  qst*B(c)
%
% mu is place-dependent (mass ratios m1/m2, m1/m5 vary along the cochlea), which
% is why it multiplies elementwise rather than being folded into Dq.
% Equivalence with the old branches verified on 11 configurations by
% macro_shadow.m -- errors 0 to 6.8e-17, including every nested/clcouple/vent
% combination and the deliberately asymmetric m=2 case.
C = cp.mc;
nd = min(size(ss,2), C.ndof);
Q  = (C.mu(:,1:nd) .* ss(:,1:nd)) * C.Dq(:,1:nd).';   % (n x nch)
Q(1,:) = Q(1,:) + qst * C.B.';        % stapes drives the basal boundary only
qx = reshape(Q.', n*C.nch, 1);        % back to place-major interleaving
end
function a=fold_p(pp,m,n,a,ss,ii,cp)
% MACRO/MICRO INTERFACE (Stage 2). Was ~45 lines of per-chamber pressure pickups
% branching on m; it is the transpose partner of xpnd_q:
%
%     s = s_internal - Df*p
%
% then each DOF's acceleration is s_k/M_k, plus a(i1) when that DOF is measured
% ABSOLUTELY (referenced to the BM). C.dref carries that per DOF -- it used to be
% keyed on chamber count, which is a micromechanics property wearing a
% macromechanics key and precisely what this refactor exists to separate.
% d2/d3 are absolute at m>=3 and RELATIVE at m<=2; the vent state is a fluid flow
% and takes no a(i1) term. Preserved exactly: changing it would break m=1==m=2.
C = cp.mc; i1 = ii(:,1);
nd = min([size(ss,2), C.ndof, size(ii,2)]);
P  = reshape(pp, C.nch, n).';                      % (n x nch), place-major
fs = ss(:,1:nd) - P * C.Df(1:nd,:).';              % pressure differences
a(i1) = fs(:,1) ./ cp.m1;
for k = 2:nd
    if (C.dref(k)), a(ii(:,k)) = a(i1) + fs(:,k) ./ dof_mass(cp,k);
    else,           a(ii(:,k)) =         fs(:,k) ./ dof_mass(cp,k);
    end
end
% pa.bmrigid = KINEMATIC clamp of the whole basilar membrane, the correct way
% to express d1 -> 0. Stiffening k1o instead raises the BM resonance and breaks
% the fixed dt=2e-6 explicit step (every stiffened run went non-finite), so that
% route measures the integrator, not the physics. A kinematic clamp adds no
% stiffness and cannot violate the time step.
% PURPOSE: species with no basilar membrane (lizards, birds, amphibians) still
% show sharp tuning and active emissions, so amplification belongs to the hair
% cell and tectorial micromechanics. A model whose micromechanics survive a
% rigid BM reproduces that; m=1/m=2 cannot, because the BM is the only partition
% the fluid acts on.
if (isfield(cp,'bmrigid') && cp.bmrigid), a(i1) = 0; else
a(i1(1)) = 0; a(i1(n)) = 0;
end
end

function M = dof_mass(cp, k)
% Mass for DOF k: 1 BM, 2 TM-RL shear, 3 OC height, 4 vent flow. The max() on
% clvm and its absence on m2/m5 reproduce the old code exactly -- do not
% "tidy" that into consistency, it would change results where a mass is tiny.
switch k
    case 2, M = cp.m2;
    case 3, M = cp.m5;
    case 4, M = max(cp.clvm, 1e-12);
    otherwise, M = cp.m1;
end
end
function [mbf,bf]=tdm_show(pa,sav,cur,ts,cp)
report_level(sav,ts);
plot_ts(pa,sav,cur,ts);
[mbf,bf]=plot_sv(pa,sav,ts,cp);
end

function pa=parm_rep(pa,pr)
% TDM parameter replacement
if (isfield(pr,'gam')), pa.gam=pr.gam;end         % CA gain
if (isfield(pr,'chl'))                            % CHL
    chl=pr.chl(1);
    pa.rma=pa.rma*chl;
    pa.kma=pa.kma*chl;
    pa.rst=pa.rst*chl;
    pa.kst=pa.kst*chl;
end
if (isfield(pr,'gampro')&&(size(pr.gampro,1)==2)) % CA-gain profile
    x=linspace(0,1,pa.n)';
    y=pr.gampro(1,:)';
    z=pr.gampro(2,:)';
    pa.gampro=interp1(z,y,x);
end
if (isfield(pr,'synpro')&&(size(pr.synpro,1)==2)) % CA-gain profile
    x=linspace(0,1,pa.n)';
    y=pr.synpro(1,:)';
    z=pr.synpro(2,:)';
    pa.synpro=interp1(z,y,x);
end
if (isfield(pr,'hbnl')), pa.hbnl=pr.hbnl; end
if (isfield(pr,'ihceq')), pa.ihceq=pr.ihceq; end
end % return

function write_data(fn,data)
[nr,nc] = size(data);
fp=fopen(fn,'wt');
fprintf(fp,'; %s\n', fn);
for i=1:nr
    for j=1:nc
        fprintf(fp,' %14.5g',data(i,j));
    end
    fprintf(fp,'\n');
end
fclose(fp);
end

% multi-tone stimulus
function [x,nn,tt]=multi_tone(stcfg,tt,nn,td)
spl_ref = 0.0002; % SPL reference pressure (rms Pa)
[nr,nf]=size(stcfg);
if (nr>1)
    cfn = 'pedcal.mat';
    if (~exist(cfn,'file'))
        error('needs eardrum-pressure calibration (%s)',cfn); end
    load(cfn,'f','pc_gn')
end
x=zeros(1,nn);
flst=stcfg(1,:);
for k=1:nf
    fr=round(flst(k)*tt)/tt; % adjust frequency for integer cycles
    if (nr==1)
        gn = 1;
    else
        lv=stcfg(2,k);
        pr=interp1(f,abs(pc_gn),fr);
        gn=(spl_ref/pr/sqrt(2))*10^((lv+40)/20);
    end
    t=gn*sin(2*pi*fr*(linspace(0,tt,nn)-td)); % compute tone
    x=x+t;                               % add tone to stimulus
end
if (nr==1)
    x=x/rms(x); % normalize rms average
end
end % return

% band-limited stimulus
function [x,nn,tt]=band_limit(stcfg,tt,nn,st,rc)
dt=tt/nn; % time-step (sec)
tr=70;    % ramp duration
trt=tr+170;
t2=trt+tr;
nn=trt/dt;
tt=nn*dt;
n2=round(t2/dt);
w=tukeywin(n2,2*tr/t2);
w=reshape(w(1:nn),1,nn); % truncate window
fr=abs(stcfg(1,1));
lv=stcfg(2,1);
bw=max(1e-6,stcfg(3,1));
svcfg=stcfg;
stcfg=[fr;lv];
[x0,nn,tt]=multi_tone(stcfg,tt,nn,0);
fr=fr*2.^(bw*[-0.50 -0.25 0 0.25 0.50]);
lv=lv*[1 1 1 1 1];
stcfg=[fr;lv];
td=-9; % time delay
if (st>0)
    [xx,nn,tt]=band_noise(svcfg,tt,nn,st,rc);
else
    [xx,nn,tt]=multi_tone(stcfg,tt,nn,td);
end
x=(xx.*w)*(rms(x0)/rms(xx));
end % return

% multi-toneburst stimulus
function [x,nn,tt]=multi_toneburst(stcfg,tt,nn)
cfn = 'pedcal.mat';
if (~exist(cfn,'file'))
    error('needs eardrum pressure calibration (%s)',cfn); end
load(cfn,'f','pc_gn')
spl_ref = 0.0002; % SPL reference pressure (rms Pa)
[nr,nf]=size(stcfg);
lvc=0.9;   % stimulus-level correction
sgn=1;     % stimulus sign
tr=10;     % ramp duration
type=nr-2; % 1=tbabr 2=fwdmsk
if (type==1)
    trt=20;
elseif(nf==1) % ECochG stimulus
    dur=abs(stcfg(3));
    sgn=sign(stcfg(3));
    trt=dur+10;
    stcfg(3)=dur;
    tr=stcfg(5); % ramp duration
else
    trt=160;
end
dt=tt/nn; % time-step (sec)
nn=round(trt/dt)+1;
tt=nn*dt;
x=zeros(1,nn);
flst=stcfg(1,:);
for k=1:nf
    fr=flst(k);
    lv=stcfg(2,k) + lvc;
    t2=stcfg(3,k);        % tone duration
    n2=round(t2/dt)+1;
    if (type==1)
        t1=0;
        w=blackman(n2)';
    else
        t1=stcfg(4,k);
        w=tukeywin(n2,2*tr/t2)';
    end
    n1=round(t1/dt)+1;
    n3=n1+n2-1;
    pr=interp1(f,abs(pc_gn),max(min(f),fr));
    gn=sgn*(spl_ref/pr/sqrt(2))*10^((lv+44)/20);
    tb=gn*cos(pi*fr*linspace(-t2,t2,n2)).*w; % compute toneburst
    if (n3>nn) % truncate tone-burst extension beyond stimulus
        tb=tb(1:(n2+nn-n3));
        n3=nn;
    end
    x(n1:n3)=x(n1:n3)+tb; % add toneburst to stimulus
end
end % return

function s=amod_tone(nt,bw,cf,sr)
nf=3;     % number of frequencies
tt=nt/sr; % total time
mf=cf*(2^bw-1)/(2^bw+1);    % modulaton frequency
flst=cf+[-mf 0 mf];         % frequency list
gain=[1 2 1];               % relative gain of crequency components
s=zeros(nt,1);              % signal
t=linspace(0,tt,nt)';       % time
for k=1:nf
    f=round(flst(k)*tt)/tt; % frequency (integer cycles)
    g=gain(k);              % gain
    c=g*sin(2*pi*f*t);     % compute tone
    s=s+c;                  % add tone to stimulus
end
%--------------------------------------------
s=s/sqrt(mean(s.^2)); % normalize to rms=1
end % return

% band-noise stimulus
function [x,nn,tt]=band_noise(stcfg,tt,nn,st,rc)
cfn = 'pedcal.mat';
if (~exist(cfn,'file'))
    error('needs eardrum-pressure calibration (%s)',cfn); end
load(cfn,'f','pc_gn')
rng('default'); % reset random number generator
if (rc>0)
    rand([1 rc]); % random-count steps from default
end
spl_ref = 0.0002; % SPL reference pressure (rms Pa)
dt=tt/nn; % time-step (sec)
tr=70;    % ramp duration
trt=tr+50;
t2=trt+tr;
nn=trt/dt;
tt=nn*dt;
lvc=-0.099; % stimulus-level correction
fr=abs(stcfg(1,1));
lv=stcfg(2,1) + lvc;
bw=abs(stcfg(3,1));
n2=round(t2/dt);
w=tukeywin(n2,2*tr/t2);
pr=interp1(f,abs(pc_gn),fr);
gn=(spl_ref/pr/sqrt(2))*10^((lv+44)/20);
if (st==1)
    noise=flat_noise(n2,bw,fr,1/dt); % compute  flat noise noise
elseif(st==2)
    noise=lono_noise(n2,bw,fr,1/dt); % compute low-noise noise
elseif(st==3)
    noise=amod_tone(n2,bw,fr,1/dt);  % compute AM tone
else
    error('undefined noise type: st=%d',st)
end
nb=gn*noise.*w; % compute noise burst
x=nb(1:nn)';     % band-noise stimulus
end % return

function ts=plot_ts(pa,sav,cur,ts,tt,smx)
persistent p1 p2 a2
if (ts.done), return; end
if (ts.tsdsp<1), return; end
if (nargin>5)      % initialize
    ts.smx = smx / 10; ts.tic = tic;
    clf; figure(ts.tsdsp)
    subplot(2,1,1); p1=plot(0,0,'-b'); xlim([0 tt]);
    title(ts.lab); xlabel('time (ms)'); ts.stim=pa.stim;
    switch(ts.dsp1), case 1, ylab='V_r'; case 2, ylab='P_e'; case 3, ylab='P_s'; case 4, ylab='WNR'; end
    ylabel(ylab); subplot(2,1,2); p2=plot(0,0,'-b'); a2=gca;
    if (pa.m), xlabel('place (mm)'); else, set(a2,'XtickLabel',[]); xlabel('stapes'); end
elseif (nargin==5) % Animation
    dspint=32; i=tt;
    if (mod(i-1,dspint)), return; end
    if (i>dspint)
        kk=1:i; xx=kk*pa.ntsw*(ts.dt);
        switch(ts.dsp1), case 1, yy=sav.vep(kk); case 2, yy=sav.ped(kk); case 3, yy=sav.pst(kk); case 4, yy=sav.wnr(kk); end
        set(p1,'XData',xx,'YData',yy)
        n = pa.n; ii = 2:(n-2); dx = 10 * pa.xl / (n - 1); x = (ii-1)' * dx;
        jj = 1:pa.ldne; xx = (jj-1)' * (10 * pa.xl / (pa.ldne - 1));
        switch(ts.dsp2), case 1, y = cur.hb(ii); case 2, y = cur.nr(ii); case 3, y = cur.ld; x = xx; case 4, y = cur.hc(ii); end
        set(p2,'XData',x,'YData',y)
        if (pa.m<1), ymxmn=2e-11; else, ymxmn=1e-9; end
        ymx=max(max(abs(y)),ymxmn);
        if (~isfinite(ymx)), ymx=ymxmn; end   % guard: Inf/NaN YLim hangs the graphics controller
        switch(ts.dsp2), case 1, ymn = -ymx; case 2, ymn = -ymx / 20; case 3, ymn = -ymx / 20; case 4, ymn = -ymx; end
        set(a2,'YLim',[ymn ymx]); drawnow;
    end
else               % Summary
    figure(ts.tsdsp); subplot(2,1,1);
    switch(ts.dsp1), case 1, yy=sav.vep; case 2, yy=sav.ped; case 3, yy=sav.pst; case 4, yy=sav.wnr; end
    nt=length(yy); xx=(1:nt)*pa.ntsw*(ts.dt); plot(xx,yy,'b')
    smx = max(ts.smx,max(abs(yy))); if (~isfinite(smx)), smx=max(ts.smx,eps); end   % same guard as line ~756
    if (ts.dsp1>3), smn=-smx/20; else, smn=-smx; end
    if (smn==smx), smn=smx-eps; end; axis([0 max(xx) [smn smx]*1.01]);
    subplot(2,1,2); kk=pa.n:-1:1; jj=pa.ldne:-1:1;
    switch(ts.dsp2), case 1, yy=sav.hb(:,kk)'; case 2, yy=sav.nr(:,kk)'; case 3, yy=sav.ld(:,jj)'; case 4, yy=sav.hc(:,kk)'; end
    xc=[0 ts.ntsf*ts.dt*10]; yc=[pa.xl*10 0]; imagesc(xc,yc,yy); xlabel('time (ms)')
    if (pa.m), ylabel('place (mm)'); else, set(gca,'YtickLabel',[]); ylabel('stapes'); end
    if (pa.m)
        if (pa.ihcv==0), hbv='displacement'; else, hbv='velocity'; end
        switch(ts.dsp2), case 1, title(sprintf('hair-bundle %s',hbv)); case 2, title('neural rate'); case 3, title('loudness spectrogram'); end
    end
    if (ts.dsp2==3)
        pld = loudness(sav.nr,pa.ldew,pa.ldne,pa.ldsc);
        ldso=mean(pld); ldph=sone2phon(ldso); fprintf('loudness = %.1f (phon)\n',ldph);
    end
end
end

function [mbf,bf]=plot_sv(pa,sav,ts,cp)
nsv=length(pa.isv); mbf=zeros(1,nsv); bf=zeros(1,nsv); % bf default for tone/early-return paths
if (ts.done || ~ts.svdsp), return; end
fig=ts.svdsp;
if (iscell(ts.stcfg)), ts.stcfg=1; end
% only a scalar-zero (click) triggers plot_click, which saves the pedcal calibration
if (~(isscalar(ts.stcfg) && ts.stcfg==0)), plot_tone(sav,fig); else, [mbf,bf]=plot_click(pa,sav,ts,cp,fig); end
end

function [mbf,bf]=plot_click(pa,sav,ts,cp,fig)
pgd=0; plot_cal(pa,sav,fig); [mbf,bf]=plot_tc(pa,sav,ts,cp,pgd,fig+1);
end

function plot_tone(sav,fig)
spl_ref = 0.0002; f=sav.f; nf=length(f); ii=1:nf; [~,~,pe]=fetch_sav(1,sav); spl=dbmag(pe/spl_ref/nf/sqrt(2));
figure(fig);clf; semilogx(f,spl(ii)); axis([0.1 20 -10 90]); xlabel('frequency (kHz)'); ylabel('P_e (dB SPL)'); title('sound level at eardrum'); drawnow
end

function plot_cal(pa,sav,fig)
ef=[0.1 0.2 0.5 1 2 5 10]; emZe=10^4*[22.8 13.3 6.5 4.8 4.1 4.0 4.9]; emZc=1e6*[8.36 6.28 9.936 10.7 15.7 22.6 19.8]; emGf=[-5.63 -2.77 11.0 16.6 11.7 0.25 -10.7];
f=sav.f; f(f<1e-9)=1e-9; nf=length(f); ii=1:nf;
[~,~,pe,ps,vr,vs,ve]=fetch_sav(1,sav);
pc_gn=pe./vr; ps_gn=ps./pe; zc=ps./(vs/pa.ast); zm=pe./ve;
f=f(ii); pc_gn=pc_gn(ii); ps_gn=ps_gn(ii); zc=zc(ii); zm=zm(ii);
if (sav.nch), save pedcal.mat f pc_gn ps_gn; end
zc=zc*10^5; zm=zm*10^5; pc_db=dbmag(pc_gn); ps_db=dbmag(ps_gn);
zc_mg=max(abs(zc),1e-9)/1000; zm_mg=max(abs(zm),1e-9)/1000;
zm_ph=angle(zm)/(2*pi); ps_ph=angle(pc_gn)/(2*pi);
figure(fig);clf; subplot(2,3,1); loglog(f,zm_mg,ef,emZe,'r'); xlim([0.1 20]); ylabel('Zm (mks k\Omega)'); title('middle-ear impedance');
subplot(2,3,4); semilogx(f,zm_ph); ylabel('phase (cyc)'); axis([0.1 20 [-1 1]*0.301]);
subplot(2,3,3); loglog(f,zc_mg,ef,emZc,'r'); xlim([0.1 20]); ylabel('Zc (mks k\Omega)'); title('cochlea impedance');
subplot(2,3,2); semilogx(f,ps_db,ef,emGf,'r'); axis([0.1 20 -40 40]); xlabel('frequency (kHz)'); ylabel('Ps / Pe (dB)'); title('middle-ear gain');
subplot(2,3,5); semilogx(f,ps_ph); ylabel('phase (cyc)'); axis([0.1 20 [0 1]*0.601]);
subplot(2,3,6);
if (pa.m)
    f=sav.f; nsv=length(pa.isv); bf=zeros(1,nsv);
    for i=1:nsv, [d1,d2]=fetch_sav(i,sav); vh=abs((d1-d2).*f); [~,ii]=max(log(vh)); bf(i)=f(ii); end
    xl=pa.xl*10; xp=xl*(pa.isv-1)/(pa.n-1); fr=logspace(-1,log10(20),100); cp=f2p(fr,0.2,0.021,0.99,xl);
    semilogx(bf,xp,'b',fr,cp,'r'); axis([0.1 20 0 xl]); xlabel('frequency (kHz)'); ylabel('place (mm)'); title('frequency-place map');
else
    semilogx(f,pc_db); axis([0.1 20 -30 30]); xlabel('frequency (kHz)'); ylabel('Pe / Vr (dB)'); title('source calibration');
end
drawnow
end

function [mbf,bf]=plot_tc(pa,sav,ts,cp,pgd,fig)
nsv=length(pa.isv); bf=zeros(1,nsv); mbf=zeros(1,nsv);
if (~pa.m), return; end; if (length(ts.stcfg)>1), return; end
f=sav.f; s=2i*pi*f; nf=length(f); 
spl=zeros(nf,nsv); phv=zeros(nf,nsv); gdv=zeros(nf,nsv); qe=zeros(1,nsv); pbf=zeros(1,nsv);
dsp_ref = 140; spl_ref = 0.0002; j1=ceil(nf*(0.1/max(f)));
for i=1:nsv
    [d1,d2,pe,~,~,~,~,nr] = fetch_sav(i,sav);
    if (pa.ihceq<1)
        vh = nr/1000; 
    else
        sv = pa.isv(i);
        % TWO-COLUMN ON PURPOSE. d1/d2 here are SAVED SPECTRA from fetch_sav, not
        % current state, so the third term would have to be the saved d3 series at
        % this place (sav.d3) transformed the same way -- fetch_sav does not
        % return it. Passing cur.d(sv+2*n) here (which I briefly did) is a scalar
        % from the FINAL timestep: wrong length and wrong quantity, and it would
        % broadcast SILENTLY because hb(:,3) is zero today.
        % So refuse rather than compute a wrong neural response.
        if (isfield(pa,'hbrl') && pa.hbrl)
            error('tdm26:hbrlIHC', ...
                ['pa.hbrl (RL-referenced bundle) is not wired into the IHC path. ' ...
                 'vh here is built from fetch_sav spectra, which do not include ' ...
                 'd3. Extend fetch_sav to return the saved d3 series before using ' ...
                 'hbrl with ihceq>=1, or set pa.ihceq<1 to take the nr/1000 path.']);
        end
        vh = hbmix(cp.hb(sv,:), d1, d2).*s;
    end
    mnvh=1e-22; vh(isnan(vh)|(abs(vh)<mnvh)) = mnvh; pe = pe / spl_ref;
    spl(:,i)=20*log10(abs(pe./vh))+pa.hbt-dsp_ref; phv(:,i)=unwrap(angle(vh./pe))/(2*pi); gdv(:,i)=delay(vh./pe,f);
    [bfi,qei,j]=find_bf(f,-spl(:,i)); bf(i)=bfi; qe(i)=qei;
    % j is NaN when find_bf rejects an edge maximum; do not index with it.
    if (isfinite(j)), mbf(i)=spl(j,i); pbf(i)=phv(j,i); else, mbf(i)=NaN; pbf(i)=NaN; end
    for j=j1:nf, if (gdv(j,i)<1e-9), phv(j:end,i)=nan; gdv(j:end,i)=nan; break; end, end
end
spl(f>2,1)=nan; spl(f>4,2)=nan;
fprintf('BF =  ');fprintf(' %5.2f', bf );fprintf('  \n');
fprintf('Qe =  ');fprintf(' %5.2f', qe );fprintf('  \n');
fprintf('pbf=-[');fprintf(' %5.2f',-pbf);fprintf('];\n');
fpk=500*2.^(-1:5); pk=pk_tgt(fpk); mpk=-pk(2,:);
figure(fig);clf; subplot(2,1,1); reset_color_index; 
semilogx(f,spl,fpk/1000,mpk,'ro',bf,mbf,'k.'); 
axis([0.1 20 0 100]); ylabel('Pe (dB SPL)');
if (pa.ihceq<1), title('NR tuning curves'); else, title('HB velocity tuning curves'); end
text(10,15,ts.lab); subplot(2,1,2); reset_color_index;
if (pgd), loglog(f,gdv); axis([0.1 20 1 30]); ylabel('delay (ms)'); text(10,20,pa.parlab);
else, semilogx(f,phv,bf,pbf,'k.'); axis([0.1 20 -9 1]); ylabel('phase (cyc)'); text(10,-7,pa.parlab); end
xlabel('frequency (kHz)'); drawnow; print('-dpng',fig,'','tc.png');
end

function pk=pk_tgt(flst)
fpk=500*2.^(-1:5); xcp= [ 0.81 0.71 0.59 0.47 0.33 0.19 0.05]; mpk=-[18.20 9.70 8.80 15.00 12.30 19.10 59.10]; ppk=-[ 3.20 3.65 3.83 3.74 3.37 2.60 1.10];
xcp=interp1(fpk,xcp,flst); mpk=interp1(fpk,mpk,flst); ppk=interp1(fpk,ppk,flst);
emZc=1e4*[8.36 6.28 9.936 10.7 15.7 22.6 19.8]; zpk=interp1(500*2.^(-1:5),emZc,flst); pk=[xcp;mpk;ppk;zpk];
end

function ts=plot_spec(~,~,ts)
nw = 4096; w = hann(nw); fs = 1 / ts.dt; y = ts.stwav;
subplot(2,1,1); t=linspace(0,ts.tt/1000,ts.nt); plot(t,y); xlim([0 max(t)]); if (ischar(ts.stcfg)), title(ts.stcfg); elseif (ts.stcfg==0), title('click'); end
subplot(2,1,2); nov=round(nw*15/16); S=spectrogram(y,w,nov,nw*2); nf=512; df=fs/nw; S=S(1:round(nf/10),:);
t=linspace(0,ts.tt/1000,size(S,2)); f=df*(0:(nf-1)); x=f2p(f,0.2,0.021,0.99,35); mxt=max(t); xlb='time (ms)';
if (mxt<1), t=t*1000; mxt=max(t); xlb='time (ms)'; end
imagesc(t,x,abs(S)); axis([0 mxt 0 35]); xlabel(xlb); ylabel('place (mm)'); title('Fourier spectrogram'); ts.done=1;
end

function report_level(sav,ts)
if (ts.done), return; end
spl_ref = 0.0002; ped=sav.ped; ped=ped - mean(ped); pmn=min(ped); pmx=max(ped); pav=rms(ped); lv1=20*log10(pav/spl_ref); lv2=20*log10((pmx-pmn)/(2*sqrt(2)*spl_ref));
fprintf('eardrum pressure: %.1f (dB SPL) %.1f (dB ppeSPL) \n',lv1,lv2);
end

% process protocol
function ts=check_protocol(stcfg,pr,st,rc)
ts.done = 0;                 % not yet done
if (iscell(stcfg))
    return;                                             % cell array ?
elseif (~ischar(stcfg))
    return;                                             % char ?
elseif (contains(stcfg,'.wav'))
    return;                                             % WAV file  ?
elseif (contains(stcfg,'tbabr'))
    ts=tbabr_protocol(pr,st,rc);                        % TBABR
elseif (contains(stcfg,'fwdmsk'))
    ts=fwdmsk_protocol(pr);                             % FWDMSK
elseif (contains(stcfg,'ecochg'))
    ts=ecochg_protocol(pr);                             % ECOCHG
elseif (contains(stcfg,'ldngrw'))
    ts=ldngrw_protocol(pr,st,rc);                       % LDNGRW
elseif (contains(stcfg,'dpoae'))
    ts=dpoae_protocol(pr);                              % DPOAE
elseif (contains(stcfg,'optimize'))
    ts=opt_protocol(pr);                                % OPTIMIZE
elseif (contains(stcfg,'tweak'))
    ts=tweak_protocol(pr);                              % TWEAK
elseif (contains(stcfg,'param'))
    ts=param_protocol(pr);                              % PARAM
elseif (contains(stcfg,'exproto'))
    ts=user_protocol(pr);                               % EXPROTO
elseif (contains(stcfg,'coupeig'))
    ts=coupeig_protocol(pr);                            % COUPEIG (Month-2 stability)
elseif (contains(stcfg,'wnr1'))
    ts=wnr1_protocol(pr);                               % single-condition WNR trace (diagnostic)
end
ts.done = 1;                 % protocol completed!
end % return

% TBABR protocol
function S=tbabr_protocol(pr,dsp1,dsp2)
if (nargin<2), dsp1=0; end   % dsp1/dsp2 (tdm26 args 3&4) gate figure generation
if (nargin<3), dsp2=0; end
fprintf('tone-burst ABR protocol...\n')
f=(0.5*2.^(0:3))';  % stimulus frequencies (kHz)
slv=[20 40 60 80];  % stimulus levels
lev=zeros(4,4);     % eardrum levels
lat=zeros(4,4);     % WNR latencies
sho=zeros(4,4);     % WNR shoulder ratio (2nd-largest peak / main; 0 = single-peaked)
oael=nan(4,4);      % tone-burst OAE latency (NaN unless pr.oae set)
oamg=nan(4,4);      % tone-burst OAE band level (dB)
abr=zeros(4,4);     % ABR wave V latency
b=12.9;c=5;d=0.413; % ABR wave V paramters
perform_calibration % calibrate, if necessary
tic;                % start stopwatch
% collect data
for j=1:4
    for k=1:4
        fr=f(j);
        lv=slv(k);
        [mlv,tpk,nch,wnr,dtms,oae,oam] = tbabr_condition(fr,lv,pr,dsp1,dsp2);
        lev(j,k)=mlv;
        lat(j,k)=tpk;
        if (isfinite(tpk)), sho(j,k)=wnr_shoulder(wnr,dtms); end   % 0 if sub-threshold
        oael(j,k)=oae; oamg(j,k)=oam;                             % tone-burst OAE (NaN unless pr.oae)
        fprintf('%5.2f %3.0f %5.1f %4.1f\n',fr,lv,mlv,tpk);
        i=lv/100;
        abr(j,k)=b*c^(-i)*fr^(-d);
    end
end
fprintf('tbabr protocol: nch=%d wall_time=%.1f s\n',nch,toc)
% analyze data
if (dsp1 || dsp2)                 % figures suppressed when both display args are 0
    figure(1);clf
    loglog(f,lat,'bo-',f,abr,'c-');
    xlabel('frequency (kHz)')
    ylabel('latency (ms)')
    axis([0.4 5 1 20])
end
write_data('tbabr.txt',[f lev lat]);
S.f=f;S.slv=slv;S.lev=lev;S.lat=lat;S.abr=abr;S.sho=sho;S.oae=oael;S.oam=oamg;
end % return

% WNR shoulder metric: ratio of the 2nd-largest peak to the main peak within a
% ~20 ms window (0 = clean single-peaked CAP, ->1 = a competing shoulder/2nd peak
% like the early second-DOF response that corrupts the Wave-I latency).
function sh=wnr_shoulder(w,dtms)
w=w(:); n=numel(w); kmax=min(n-1,max(2,round(20/dtms))); ii=(2:kmax)';
lm=ii(w(ii)>w(ii-1) & w(ii)>=w(ii+1)); p=w(lm); p=p(p>0);
if (numel(p)<2), sh=0; return; end
p=sort(p,'descend'); sh=p(2)/p(1);
end

% forward-masking  protocol
function S=fwdmsk_protocol(pr)
if (pr), S=fwdmsk_subprotocol; return; end
fprintf('forward-masking protocol...\n')
mlv=[40 60 80];    % masker level (dB SPL)
plv=20;            % probe level (dB SPL)
dly=20;            % probe delay (ms)
trt=160;           % total response time (ms)
dt=0.02;           % response time step (ms)
prb=zeros(length(mlv),length(dly));
col='rgbk';
perform_calibration % calibrate, if necessary
tic;                % start stopwatch
figure(2);clf       % clear figure 2
fprintf('; mlv plv dly prbmax\n');
for j=1:length(mlv)
    for k=1:length(dly)
        wnr=fwdmsk_condition(mlv(j),plv,dly(k));
        nt=length(wnr);
        t=linspace(0,trt,nt)';
        fn=sprintf('fwdmsk%d%d.txt',j,k);
        write_data(fn,[t wnr]);
        figure(2); plot(t,wnr,col(j)); hold on
        iprb=round((100+dly(k))/dt):nt;
        prbmax=max(wnr(iprb));
        fprintf('  %3.0f %3.0f %3.0f %8.3g\n',mlv(j),plv,dly(k),prbmax);
        prb(j,k)=prbmax;
    end
end
% unmasked probe
wnr=fwdmsk_condition(-90,plv,20);
write_data('fwdmsk00.txt',[t wnr]);
fprintf('  %3.0f %3.0f %3.0f %8.3g\n',-90,plv,20,prbmax);
fprintf('fwdmsk protocol: elapsed time = %.3f (sec)\n',toc)
% plot results
figure(2); plot(t,wnr,col(4)); hold off
xlabel('time (ms)')
ylabel('WNR')
legend('40','60','80','-90')
iprb=round((100+20)/dt):round(trt/dt);
prbmax=max(wnr(iprb));
S.mlv=mlv; S.plv=plv; S.dly=dly; S.prb=prb; S.prbmax=prbmax;
end % return

% forward-masking  subprotocol
function S=fwdmsk_subprotocol
fprintf('forward-masking subprotocol...\n')
mlv=[40 60 80]';    % masker level (dB SPL)
dly=[10 20 40]';    % probe delay (ms)
tlv=20;            % probe level target (dB SPL)
amt0=[15.59    11.92     8.25   % 40
      25.10    19.19    13.28   % 60
      34.61    26.46    18.31]; % 80
dt=0.02;           % response time step (ms)
prb=zeros(length(mlv),length(dly));
write_data('amtmsk0.txt',[mlv amt0]); % measured amount of masking
perform_calibration % calibrate, if necessary
tic;                % start stopwatch
fprintf('; mlv dly plv\n');
% unmasked probe
wnr=fwdmsk_condition(-90,tlv,dly(2));
nt=length(wnr);
iprb=round((100+dly(2))/dt):nt;
wnrref=max(wnr(iprb));
maxdev=6;
for j=1:length(mlv)
    for k=1:length(dly)
        plv=tlv+amt0(j,k)+[-1 0 1]*maxdev;
        wnrmax=zeros(size(plv));
        for i=1:length(plv)
            wnr=fwdmsk_condition(mlv(j),plv(i),dly(k));
            wnrmax(i)=max(wnr(iprb));
        end
        if (min(wnrmax)>=wnrref),     prblev=min(plv);
        elseif (max(wnrmax)<=wnrref), prblev=max(plv);
        end
        for i=2:length(plv)
            if ((wnrmax(i-1)<=wnrref)&&(wnrmax(i)>wnrref))
                a=(wnrmax(i)-wnrref)/(wnrmax(i)-wnrmax(i-1));
                prblev=a*plv(i-1)+(1-a)*plv(i);
                break;
            end
        end
        prb(j,k)=prblev;
        fprintf('  %3.0f %3.0f %3.0f\n',mlv(j),dly(k),prb(j,k));
    end
end
fprintf('fwdmsk subprotocol: elapsed time = %.3f (sec)\n',toc)
% plot results
amt=prb-tlv; % model amount of masking
figure(2); plot(mlv,amt0); hold off
xlabel('masker level (dB SPL)')
ylabel('amount of masking (dB)')
legend('10','20','40')
maxdev=max(abs(amt(:)-amt0(:)));
fprintf('fwdmsk maxdev=%.1f\n',maxdev)
write_data('amtmsk.txt',[mlv amt]);
S.mlv=mlv; S.plv=plv; S.dly=dly; S.tlv=tlv; S.prb=prb;
end % return

% ECochG  protocol
function S=ecochg_protocol(pr)
fprintf('ECochG protocol...\n')
nfr=length(pr.fr);  % number of frequencies
npo=length(pr.xe);  % number of conditions
am=zeros(nfr,npo);
gd=zeros(nfr,npo);
ph=zeros(nfr,npo);
if (~isfield(pr,'fr')) , pr.fr=0.5;  end % tone frequency (kHz)
if (~isfield(pr,'lv')) , pr.lv=60;   end % tone level (dB SPL)
if (~isfield(pr,'td')) , pr.td=20;   end % stimulus duration (ms)
if (~isfield(pr,'tr')) , pr.tr=4;    end % ramp duration (ms)
if (~isfield(pr,'ts')) , pr.ts=0;    end % stimulus start (ms)
if (~isfield(pr,'tp')) , pr.tp=40;   end % plotted time (ms)
if (~isfield(pr,'sd')) , pr.sd=0.25; end % synaptic delay (ms)
if (~isfield(pr,'xe')) , pr.xe=25;   end % electrode positions (mm)
if (~isfield(pr,'hcc')), pr.hcc=2e4; end % hair-cell current multiplier
if (~isfield(pr,'nrc')), pr.nrc=0;   end
if (~isfield(pr,'gam')), pr.gam=1;   end
if (~isfield(pr,'wd')) , pr.wd=1;    end
if (~isfield(pr,'hca')), pr.hca=0.999; end
if (~isfield(pr,'wvc')), pr.wvc = 1.5; end % Standard decay
if (~isfield(pr,'sumlim')), pr.sumlim=0; end
if (~isfield(pr,'diflim')), pr.diflim=0; end
perform_calibration % calibrate, if necessary
prf=pr;
lg=cell(1,nfr);
for kfr=1:nfr
    prf.fr=pr.fr(kfr);
    lg{kfr}=sprintf('%.2f',prf.fr);
    [ecochg0,prf]=ecochg_condition(prf,0);
    [ecochg1,prf]=ecochg_condition(prf,1);
    fprintf('ecochg stimulus: frequency=%.2f kHz\n',prf.fr)
    % compute waveforms
    rs=ecochg0+ecochg1;
    rd=ecochg0-ecochg1;
    nt=size(rd,1);
    nf=1+nt/2;
    Rs=ffa(rs);
    Rd=ffa(rd);
    aRs=abs(Rs)*2/nt;
    aRd=abs(Rd)*2/nt;
    nx=length(prf.zvc);
    tt=prf.tt;
    t=linspace(0,tt,nt)';
    f=linspace(0,(nf-1)/tt,nf);
    if (pr.wd) % write data ???
        write_data('ecochg.txt',[t rs rd]);
    end
    % save first-electrode waveforms
    if (kfr==1)
        rds=zeros(nt,nfr);
        rss=zeros(nt,nfr);
    end
    rds(:,kfr)=rd(:,1);
    rss(:,kfr)=rs(:,1);
    % plot results
    figure(2);clf       % clear figure
    subplot(2,2,1)
    plot(t,rd)
    xlim([0 pr.tp])
    if (length(pr.sumlim)>1), ylim(pr.sumlim); end
    ylabel('difference amplitude (\muV)')
    title('waveform')
    subplot(2,2,3)
    plot(t,rs)
    xlim([0 pr.tp])
    if (length(pr.diflim)>1), ylim(pr.diflim); end
    xlabel('time (ms)')
    ylabel('summed amplitude (\muV)')
    subplot(2,2,2)
    semilogx(f,aRd)
    xlim([0.1 max(4, max(pr.fr)*1.5)])
    ylabel('difference magnitude (\muV)')
    title('spectrum')
    subplot(2,2,4)
    semilogx(f,aRs)
    xlim([0.1 4])
    xlabel('frequency (kHz)')
    ylabel('summed magnitude (\muV)')
    drawnow
    for kpo=1:npo
        wt=rd(:,kpo).^2;
        [mx,ix]=max(aRd(:,kpo));
        am(kfr,kpo)=mx;
        ph(kfr,kpo)=angle(Rd(ix,kpo));
        gd(kfr,kpo)=(t'*wt)/sum(wt);
        if (pr.wd) % write data ???
            fprintf('position(%2d)=%.1f ',kpo,pr.xe(kpo))
            fprintf('amplitude=%4.1f ',am(kfr,kpo))
            fprintf('delay=%4.1f ',gd(kfr,kpo))
            fprintf('shift=%4.1f ms\n',gd(kfr,kpo)-gd(kfr,1))
        end
    end
    figure(3);clf % clear figure
    xl=prf.xl;
    x=linspace(0,xl,nx);
    subplot(2,1,1)
    plot(x,prf.zvc)
    xlim([0 xl])
    xlabel('cochlear place (mm)')
    ylabel('VC impedance')
    title('electrode location')
    subplot(2,1,2)
    plot(t,rd); hold on
    reset_color_index;
    plot(gd(:), zeros(size(gd(:))), 'o', 'LineWidth', 2)
    xlim([0 pr.tp])
    if (length(pr.diflim)>1), ylim(pr.diflim); end
    ylabel('ECochG difference')
    xlabel('time (ms)')
    drawnow
end
S.xe=pr.xe;
S.fr=pr.fr;
S.td=pr.td;
S.tt=tt;
S.rd=rds;
S.rs=rss;
S.am=am;
S.gd=gd;
S.ph=ph;
S.lg=lg;
if (nfr>2)
    figure(4);clf % clear figure
    ii=am<0.6;    % don't plot amplitudes < 0.6 microvolt
    am(ii)=nan;
    gd(ii)=nan;
    xe=pr.xe;
    xx=[min(xe)-0.5 max(xe)+0.5];
    subplot(2,1,1)
    plot(xe',am','-o')
    xlim(xx)
    ylabel('amplitude (\muV)')
    title('ECochG difference')
    subplot(2,1,2)
    plot(xe',gd','-o')
    xlim(xx)
    ylabel('delay (ms)')
    xlabel('electrode location (mm)')
    drawnow
end
end % return

% LDNGRW - loudness-growth protocol
function S=ldngrw_protocol(pr,st,rc)
fprintf('loudness growth protocol: ')
slv=20*(0:5)';                % stimulus level
sbw=[0.0001 [1 2 3 4 8]/8]';  % stimulus bandwidth (octaves)
sfr=1;                        % stimulus frequencies (kHz)
if (pr==1)
    fprintf('loudness summation ')
    fprintf('(st=%d rc=%d)\n',st,rc)
elseif (pr==2)
    fprintf('ensemble width ')
    fprintf('(st=%d rc=%d)\n',st,rc)
    S=ldnenw_subprotocol(st,rc);
    return;
elseif (pr==3)
    fprintf('number of ensembles ')
    fprintf('(st=%d rc=%d)\n',st,rc)
    S=ldnnen_subprotocol(st,rc);
    return;
elseif (pr==4)
    fprintf('loudness summation @ 60 dB SPL ')
    fprintf('(st=%d rc=%d)\n',st,rc)
    slv=60;              % stimulus level
else
    fprintf('%.1f-kHz tone\n',sfr)
    sbw=0;               % single-tone only
end
bwk=sfr*(2.^(sbw/2)-2.^(-sbw/2)); % stimulus bandwidth (kHz)
bwk(bwk<1e-6)=1e-6;
%--------------------------------------------
nlv=length(slv);
nbw=length(sbw);
ldn=zeros(nlv,nbw);
lcu=zeros(nlv,nbw);
lph=zeros(nlv,nbw);
perform_calibration % calibrate, if necessary
pa=modpar26(1);
fprintf('ensemble: ldew=%.1f; ldne=%d\n',pa.ldew,pa.ldne);
fprintf('  fr   lv   bw   ldph\n');
tic;                % start stopwatch
% collect data
for j=1:nbw
    bw=sbw(j);
    dn=sprintf('bw%03d',round(abs(bw),2)*100);
    savnr=exist(dn,'dir')&&(st==1);
    savvep=exist('vep','dir');
    for k=1:nlv
        lv=slv(k);
        sav = ldngrw_condition(sfr,lv,bw,st,rc);
        ldso = sav.ldso;
        ldcu = sone2cu(ldso);
        ldph = sone2phon(ldso);
        ldn(k,j) = ldso;
        lcu(k,j) = ldcu;
        lph(k,j) = ldph;
        fprintf('%5.2f %3.0f %5.2f %5.2f\n',sfr,lv,bw,ldph);
        if (savnr)
            fn=fulllfile(dn,sprintf('ldngrw%03d.mat',lv));
            nr=sav.nr;
            fr=sfr;
            save(fn,'nr','fr','bw','bwk');
        end
        if (savvep)
            bwpo=round(bw*100); % bandwidth (% oct)
            fn=fullfile('vep',sprintf('vep_%03d_%03d.mat',lv,bwpo));
            vep=sav.vep;
            tt=sav.tt;
            save(fn,'vep','tt');
        end
    end
end
fprintf('ldngrw protocol: elapsed time = %.3f (sec)\n',toc)
% CLS data - Rasetshwane et al. 2015
lvd=[0.17 20.14 40.10 60.07 80.03 100.00];
ldd=[0.96  6.64 11.69 18.00 25.19  47.15];
ldp=cu2phon(ldd);
write_data('ldnclsnh.txt',[lvd' ldd' ldp']);
% Fletcher & Munson
[son,bet]=fm33;
% plot model data
figure(1); clf
if (pr==0)
    subplot(1,2,1)
    semilogy(slv,ldn,bet,son,'k:')
    xlabel('level (dB SPL)')
    ylabel('loudness (sone)')
    legend('tdm26','Fletcher','Location','SouthEast')
    axis([-5 105 0.001 100])
    subplot(1,2,2)
    plot(slv,lph,lvd,ldp,'k:')
    xlabel('level (dB SPL)')
    ylabel('loudness (phon)')
    axis([-5 105 -5 105])
    legend('tdm26','Rasetshwane','Location','NorthWest')
    drawnow
    %----
    write_data('ldngrwlv.txt',[slv ldn lcu lph]);
    save ldncal slv ldn lcu lph
    sbw=0;
elseif (pr==4) % @ 60 dB SPL
    plot(slv,lph,'+',lvd,ldp,'k:')
    xlabel('level (dB SPL)')
    ylabel('loudness (phon)')
    axis([-10 110 -10 110])
    legend('0','1/4','1','Location','Best')
    drawnow
    qolr=1-(lph(2)/lph(1));
    fprintf('qolr=%.1f%%\n',qolr*100);
    %----
    elv=eqvlvl(lph);
    write_data('ldngrwbw.txt',[bwk lcu' elv' lph']);
else                     % loudness summation for level & bandwidth
    subplot(1,2,1)
    plot(slv,lph,lvd,ldp,'k:')
    xlabel('level (dB SPL)')
    ylabel('loudness (phon)')
    axis([-10 110 -10 110])
    legend('0','1/4','1','Location','Best')
    subplot(1,2,2)
    plot(sbw,lph')
    xlabel('bandwith (oct)')
    ylabel('loudness (phon)')
    axis([-0.05 1.05 -10 110])
    legend('0','20','40','60','80','100','Location','Best')
    drawnow
    %----
    write_data('ldnsumlv.txt',[slv lph]);
    write_data('ldnsumbw.txt',[bwk lph']);
end
S.slv=slv; S.ldn=ldn; S.lcu=lcu; S.lph=lph; S.sbw=sbw;
end % return

function elv=eqvlvl(ldph)
cfn='ldncal.mat';
if (~exist(cfn,'file'))
    error('needs loudness-growth calibration (%s)',cfn); end
load(cfn,'slv','lph') % use previous loudness function for calibration
nc=length(ldph);
elv=zeros(size(ldph)); % equivalent level
for k=1:nc
    elv(k) = interp1(lph,slv,ldph(k)); % convert CU to dB SPL
end
end % return

% LDNENW protocol
function S=ldnenw_subprotocol(st,rc)
cfn='ldncal.mat';
if (~exist(cfn,'file'))
    error('needs loudness-growth calibration (%s)',cfn); end
load(cfn,'slv','lcu') % use previous loudness function for calibration
clv=slv;              % calibration levels
sbw=[0.001 [1 2 3 4 8]/8]';  % stimulus bandwidth (octaves)
enw=[1 2 3]';         % ensemble width (mm)
sfr=1;                % stimulus frequency (kHz)
slv=60;               % stimulus level
bwk=sfr*(2.^(sbw/2)-2.^(-sbw/2)); % stimulus bandwidth (Hz)
bwk(bwk<0.005)=0.005;
ne=0;
nr=length(sbw);
nc=length(enw);
elv=zeros(nr,nc);    % equivalent level
% collect data
fprintf('  sfr  slv  bw    ew   elv\n');
for j=1:nr
    for k=1:nc
        bw=sbw(j);
        ew=enw(k);
        sav = ldngrw_condition(sfr,slv,bw,st,rc,ew,ne);
        elv(j,k) = interp1(lcu,clv,sav.ldcu); % convert CU to dB SPL
        elv(j,k) = elv(j,k) - elv(1,k) + slv; % match tone level
        fprintf('%5.2f %3.0f %5.2f %5.2f %5.1f\n',sfr,slv,bw,ew,elv(j,k));
    end
end
% save data
write_data('ldnenw.txt',[bwk elv]);
S.slv=slv; S.elv=elv; S.sbw=sbw; S.bwk=bwk; S.enw=enw;
end % return

% LDNNEN protocol
function S=ldnnen_subprotocol(st,rc)
cfn='ldncal.mat';
if (~exist(cfn,'file'))
    error('needs loudness-growth calibration (%s)',cfn); end
load(cfn,'slv','lcu') % use previous loudness function for calibration
clv=slv;              % calibration levels
sbw=[0.001 [1 2 3 4 8]/8]';  % stimulus bandwidth (octaves)
nen=[35 70 140]';     % number of ensembles
sfr=1;                % stimulus frequency (kHz)
%slv=55.7;             % stimulus level
slv=60;               % stimulus level
bwk=sfr*(2.^(sbw/2)-2.^(-sbw/2)); % stimulus bandwidth (kHz)
bwk(bwk<0.005)=0.005;
ew=0;
nr=length(sbw);
nc=length(nen);
elv=zeros(nr,nc);    % equivalent level
% collect data
for j=1:nr
    for k=1:nc
        bw=sbw(j);
        ne=nen(k);
        sav = ldngrw_condition(sfr,slv,bw,st,rc,ew,ne);
        elv(j,k) = interp1(lcu,clv,sav.ldcu); % convert CU to dB SPL
        elv(j,k) = elv(j,k) - elv(1,k) + slv; % match tone level
        fprintf('%5.2f %3.0f %5.2f %5.0f %5.1f\n',sfr,slv,bw,ne,elv(j,k));
    end
end
% save data
write_data('ldnnen.txt',[bwk elv]);
S.slv=slv; S.elv=elv; S.sbw=sbw; S.bwk=bwk; S.nen=nen;
end % return

% DPOAE protocol
function S=dpoae_protocol(~)
fprintf('DPOAE protocol: not yet implemented\n');
S.sav=0;
end % return

% PARAM protocol
function S=param_protocol(rpa)
% default to 1D click response
stcfg = 0;
nch = 1;
dsp1 = 1;
dsp2 = 2;
ts.done = 0;
% compute time-domain model response ------------
[~,sav,cur,ts] = tdm_init(stcfg,nch,dsp1,dsp2,ts); % initialize model
pa = rpa;                                          % replace parameters
[sav,cur,ts,cp] = tdm_step(pa,sav,cur,ts);         % step through time
       [mbf,bf] = tdm_show(pa,sav,cur,ts,cp);      % plot saved variables
%-------------------------------
[~,~,pe,ps,~,~,ve]=fetch_sav(1,sav);
sav.gme=ps./pe; % middle-ear gain
sav.zme=pe./ve; % acoustic impedance
%-------------------------------
sav.pa = pa;   % save model parameters
sav.ts = ts;   % save time-step variables
sav.cp = cp;   % save cochlear-partition variables
sav.bf = bf;   % save best frequencies
sav.mbf = mbf; % save threshold levels
%-------------------------------
S = sav;     % return structure of saved variables
end % return

function phon=cu2phon(CU)
% reference phon for NH1k listeners from Rasetshwane et al. (2015)
NH1k=[0 18.309 33.668 48.953 63.718 75.942 82.27 87.093 91.878 96.643 102];
CU(CU<00)=00;
CU(CU>50)=50;
phon=zeros(size(CU));
for k=1:length(CU)
    phon(k)=interp1(0:5:50,NH1k,CU(k));
end
phon=round(phon,3);
end

% sone-to-CU conversion (Jesteadt et al. 2017)
function cu=sone2cu(sone)
cu=2.6253*log10(sone+0.0887).^3+0.7799*log10(sone+0.0887).^2+8.0856*log10(sone+0.0887)+13.4493;
end

% sone-to-phon conversion
function phon=sone2phon(sone)
[son,bet]=fm33;
sone=max(son(1),min(sone,son(end)));
phon=zeros(size(sone));
for k=1:length(sone)
    phon(k)=interp1(son,bet,sone(k));
end
phon=round(phon,3);
end % return

% Fletcher & Muson 1933 (sone)
function [son,bet]=fm33
bet=-10:129;
flu=[0.015     0.025      0.04      0.06      0.09      0.14      0.22      0.32 ...
      0.45       0.7         1       1.4       1.9      2.51       3.4      4.43 ...
       5.7      7.08         9      11.2      13.9      17.2      21.4      26.6 ...
      32.6      39.3      47.5      57.5      69.5      82.5      97.5       113 ...
       131       151       173       197       222       252       287       324 ...
       360       405       455       505       555       615       675       740 ...
       810       890       975      1060      1155      1250      1360      1500 ...
      1640      1780      1920      2070      2200      2350      2510      2680 ...
      2880      3080      3310      3560      3820      4070      4350      4640 ...
      4950      5250      5560      5870      6240      6620      7020      7440 ...
      7950      8510      9130      9850  1.06e+04  1.14e+04  1.24e+04  1.35e+04 ...
  1.46e+04  1.58e+04  1.71e+04  1.84e+04  1.98e+04  2.14e+04  2.31e+04   2.5e+04 ...
  2.72e+04  2.96e+04  3.22e+04   3.5e+04   3.8e+04  4.15e+04   4.5e+04   4.9e+04 ...
   5.3e+04   5.7e+04   6.2e+04  6.75e+04   7.4e+04   8.1e+04   8.8e+04   9.7e+04 ...
  1.06e+05  1.16e+05  1.26e+05  1.38e+05   1.5e+05  1.64e+05   1.8e+05  1.97e+05 ...
  2.15e+05  2.35e+05   2.6e+05  2.88e+05  3.16e+05  3.46e+05   3.8e+05  4.18e+05 ...
   4.6e+05  5.06e+05  5.56e+05  6.09e+05  6.68e+05  7.32e+05     8e+05  8.75e+05 ...
  9.56e+05 1.047e+06  1.15e+06 1.266e+06];
son=flu/flu(51); % convert FLU to sone
end % return

% execute stimulus calibration
function sav=perform_calibration
if (exist('pedcal.mat','file')), return; end
fprintf('stimulus calibration...\n')
ts.done = 0;
ts.tsdsp = 0;   % time-step figure number
ts.svdsp = 0;   % saved-variables figure number
[pa,sav,cur,ts]=tdm_init(0,1,1,1,ts);    % initialize model
pa.ihceq=0;     % no synapse
pa.hbnl=0;      % no nonlinearity
[sav,cur,ts,cp] = tdm_step(pa,sav,cur,ts);          % step through time
                  tdm_show(pa,sav,cur,ts,cp);       % plot saved variables
end % return

% execute TBABR stimulus condition
% Single-condition WNR trace for latency-detector diagnostics.
function S=wnr1_protocol(pr)
[S.mlv,S.tpk,S.nch,S.wnr,S.dtms,S.oae,S.oam,S.od,S.dgn]=tbabr_condition(pr.fr,pr.lv,pr.pa,0,0);
S.done=1;
end

function [mlv,tpk,nch,wnr,dtms,oae,oam,od,dgn]=tbabr_condition(fr,lv,pr,dsp1,dsp2)
if (nargin<4), dsp1=0; end       % dsp1/dsp2=0 suppresses the per-step animation
if (nargin<5), dsp2=0; end
if (isstruct(pr))
    pa=pr; nch=pa.m;
else
    if (isscalar(pr)), nch=max(1,pr); else, nch=1; end
    pa=modpar26(nch);
end
ts.done = 0;
spl_ref = 0.0002; % SPL reference pressure (rms Pa)
td=4/sqrt(fr); % tone-burst duration (sec)
pa0=pa;                                          % pre-init base for optional smooth OAE reference
oae_on=isfield(pa,'oae') && ~isempty(pa.oae) && pa.oae;
if (oae_on && ~isfield(pa,'rough_amp')), pa.rough_amp=3e-2; end  % OAE needs roughness above numerical floor (>=1e-2)
if (oae_on && isfield(pa,'rough_percf') && pa.rough_percf), pa.rough_fc=fr; end  % per-CF roughness window (diagnostic probe)
[pa,sav,cur,ts]=tdm_init([fr;lv;td],nch,dsp1,dsp2,ts,pa); % initialize model
pa.ihceq=4;     % Neely synapse model
pa.hbnl=1;      % nonlinear OHC transduction
            sav=tdm_step(pa,sav,cur,ts);       % step through time
ped=sav.ped;
ped=ped - mean(ped);      % remove DC
pmn=min(ped);
pmx=max(ped);
mlv=20*log10((pmx-pmn)/(2*sqrt(2)*spl_ref));
dtms=ts.dt*ts.nw;
tpk=wnr_latency(sav.wnr,dtms);          % onset-peak latency (robust to late artifacts)
wnr=sav.wnr;                            % expose the WNR trace for diagnostics
% diagnostic: peak partition displacements at the saved places. d2 is what the
% OHC compression law tests against hbmx (hbt=max(|d2|/hbmx,1)), so if
% dgn.d2mx << pa.hbmx the gain never compresses and the model runs LINEAR.
dgn.d1mx=max(abs(sav.d1(:))); dgn.d2mx=max(abs(sav.d2(:)));
dgn.d3mx=max(abs(sav.d3(:))); dgn.d3fin=all(isfinite(sav.d3(:)));
dgn.hbmx=pa.hbmx; dgn.ratio=dgn.d2mx/pa.hbmx;
% place of the dominant neural response: if this jumps basalward abruptly with
% level, the "compression cliff" is a PLACE SHIFT, not a gain effect.
enr=sum(abs(sav.nr),1); [~,ipk]=max(enr);
dgn.pkplace=ipk; dgn.pkfrac=ipk/pa.n;
w=enr/max(sum(enr),eps); dgn.pkcen=sum((1:pa.n).*w)/pa.n;   % normalized place centroid
dgn.ohcNaN=sum(~isfinite(sav.ohcp)); dgn.ohcN=numel(sav.ohcp);
dgn.ohcP=mean(sav.ohcp(isfinite(sav.ohcp)));                 % >0 = energy injected
dgn.ohcBM=mean(sav.ohcbm(isfinite(sav.ohcbm)));              % >0 = BM amplified
dgn.ohcBMW=sum(sav.ohcbm(isfinite(sav.ohcbm)))*dtms/1000;    % integrated BM work
dgn.ohcW=sum(sav.ohcp(isfinite(sav.ohcp)))*dtms/1000;
oae=NaN; oam=NaN; od=[];                % tone-burst OAE latency/level/diag (NaN/[] unless pa.oae)
if (oae_on)                             % emission = rough - smooth (stimulus & passive m.e. cancel)
    pas=pa0; if (isfield(pas,'rough_amp')), pas=rmfield(pas,'rough_amp'); end
    tss.done=0;
    [pas,savs,curs,tss]=tdm_init([fr;lv;td],nch,dsp1,dsp2,tss,pas);
    pas.ihceq=4; pas.hbnl=1;
    savs=tdm_step(pas,savs,curs,tss);
    peds=savs.ped-mean(savs.ped);
    [oae,oam,od]=tboae_latency(ped,peds,fr,dtms);
end
end % return

% onset-peak latency of the whole-nerve response
function tpk=wnr_latency(w,dtms)
% Latency = time of the LARGEST whole-nerve-response peak (the Wave-I analog)
% within a physiological window. The WNR is a modeled compound action potential
% (~ABR Wave I); the measured data latency is the Wave-V peak minus a fixed 5 ms
% neural (brainstem I->V) delay = cochlear travel time, which the nerve-level WNR
% already excludes, so NO offset is applied to the model here. Taking the largest
% LOCAL maximum (not the first threshold crossing) matches "peak of the response"
% and is robust to (a) a smaller EARLY shoulder -- the failure mode of a first-
% local-max rule, which a prior fit EXPLOITED by growing a spurious early bump to
% report an artificially short latency -- and (b) a late slow numerical ramp in
% near-critical configs, which has no interior local max. NaN when no peak exceeds
% 3x the pre-onset baseline (sub-threshold). The model WNR carries no cochlear
% microphonic, so a genuine WNR should be single-peaked (no physiological shoulder).
w=w(:); n=numel(w);
b0=mean(w(1:max(1,round(0.02*n))));      % pre-onset baseline (first 2%)
kmax=min(n-1, max(2,round(20/dtms)));    % physiological window (~20 ms; excl. late ramp)
ii=(2:kmax)';
lm=ii(w(ii)>w(ii-1) & w(ii)>=w(ii+1));   % interior local maxima
if (isempty(lm)), tpk=NaN; return; end
[pk,jj]=max(w(lm)); ix=lm(jj);           % the LARGEST peak = Wave-I peak
if (~isfinite(pk) || pk<3*b0), tpk=NaN; return; end
tpk=(ix-1)*dtms;
end % return

% COUPLED linearized stability of the partition+fluid operator (Month-2 scoping).
% Reuses the EXACT accel() (fluid solve included) as a matrix-free linear operator
% A on state [d;v]: d/dt[d;v] = [v; accel(d,v)] with stm=0 and hbnl=0 (small-signal).
% Rightmost eig(A): Re>0 => the COUPLED model is physically unstable (an implicit
% integrator cannot help); Re<0 with the explicit scheme still diverging => the
% instability is numerical (implicit/finer-dt would help).
function S=coupeig_protocol(pr)
if (~isstruct(pr)), pr=struct('nch',pr); end
nch=1; boost=0; bfrac=0.35;
if (isfield(pr,'nch')),       nch=pr.nch; end
if (isfield(pr,'boost')),     boost=pr.boost; end
if (isfield(pr,'basalfrac')), bfrac=pr.basalfrac; end
if (isfield(pr,'pa')), pa=pr.pa; nch=pa.m; else, pa=modpar26(nch); end % pr.pa: test a modified model
if (~isfield(pa,'gampro')||isempty(pa.gampro)), pa.gampro=ones(pa.n,1); end
xf=((0:pa.n-1)')/(pa.n-1);
pa.gampro = pa.gampro .* (1 + boost*0.5.*(1+cos(pi*min(xf/bfrac,1))));   % basal boost (matches gampro sweep)
pa.hbnl=0;                                   % linear / small-signal
dof = resolve_dof(pa);   % SAME function as tdm_init -- the two used to be
                         % duplicated and had to be kept in sync by hand.
pa.dof=dof; pa.ncp=pa.n*dof;
if (pa.nmev<1), pa.nmev=1; end
nsv=pa.ncp+pa.nmev; N=2*nsv;
me=midear(pa); [cp,me]=macro26(pa,me);
% Build the coupled operator densely (eigs 'largestreal' does not converge for
% this interior-eigenvalue spectrum), one accel/fluid-solve per column.
tb=tic; A=zeros(N,N);
for j=1:N
    e=zeros(N,1); e(j)=1;
    A(:,j)=coupeig_apply(e,pa,cp,me,nsv);
end
wantvec = isfield(pr,'eigvec') && pr.eigvec;   % opt-in: also return eigenvectors
if (wantvec)
    te=tic; [Vm,Dm]=eig(A); lam=diag(Dm); teig=toc(te); tbld=toc(tb)-teig;
else
    te=tic; lam=eig(A);      teig=toc(te); tbld=toc(tb)-teig;
end
[~,ord]=sort(real(lam),'descend'); lam=lam(ord);
if (wantvec), Vm=Vm(:,ord); end
S.lam=lam; S.nch=nch; S.boost=boost; S.done=1;
% spurious perturbed-zero modes come from constrained boundary DOFs (clamped BM
% ends + zero-mass helicotrema); flag PHYSICAL modes as oscillatory (f>1 kHz).
osc = abs(imag(lam))/2/pi > 1000;
S.maxRe = real(lam(1)); S.maxRe_osc = max(real(lam(osc)));
if (wantvec)
    % Decompose the MOST-UNSTABLE eigenvector (max real part) by DOF and place,
    % to identify WHICH coordinate the divergence occupies. State = [d; v]; the
    % displacement block splits into BM(1:n), shear(n+1:2n), OC-height(2n+1:3n),
    % then middle-ear evars. eig returns unit-norm vectors, so per-DOF norms
    % partition participation. Peak place x/L = ix/n (higher = more apical).
    n=pa.n; dv=pa.dof; u=Vm(:,1); ud=u(1:dv*n);
    S.uvec=u; S.ulam=lam(1);
    S.uDOFfrac=zeros(1,dv); S.uDOFpkx=zeros(1,dv);
    tot=norm(ud)+eps;
    for kk=1:dv
        seg=ud((kk-1)*n+(1:n)); S.uDOFfrac(kk)=norm(seg)/tot;
        [~,ix]=max(abs(seg));   S.uDOFpkx(kk)=ix/n;
    end
    S.uMEfrac=norm(u(dv*n+1:nsv))/tot;   % middle-ear participation
end
fprintf('coupeig nch=%d boost=%.3f  (build %.1fs eig %.1fs N=%d)\n', nch, boost, tbld, teig, N);
fprintf('  top modes [Re/s, f kHz]: ');
for t=1:8, fprintf('[%+.3g, %.1f] ', real(lam(t)), abs(imag(lam(t)))/2/pi/1000); end
fprintf('\n  maxRe(all)=%+.4g  maxRe(osc>1kHz)=%+.4g  => %s\n', ...
    S.maxRe, S.maxRe_osc, ternstr(S.maxRe_osc>0,'RHP PHYSICALLY UNSTABLE','LHP stable'));
end

function s=ternstr(c,a,b), if c, s=a; else, s=b; end, end

function Ay=coupeig_apply(y,pa,cp,me,nsv)
st.d=y(1:nsv); st.v=y(nsv+1:2*nsv); st.stm=0;
[~,a]=accel(pa,cp,me,st,st);
a(~isfinite(a))=0;   % constrain the zero-mass helicotrema DOF (m2(n)=0 -> 0/0); fully decoupled
Ay=[st.v; a];
end

% execute forward-masking stimulus condition
function wnr=fwdmsk_condition(mlv,plv,dly)
ts.done = 0;
mfr=1;       % masker frequency (kHz)
mt1=100;     % masker duration (msec)
mt2=0;       % masker start    (msec)
pfr=1;       % probe frequency (kHz)
pt1=20;      % probe duration (msec)
pt2=mt1+dly; % start of probe (msec)
stcfg=[mfr pfr;mlv plv;mt1 pt1;mt2 pt2];
[pa,sav,cur,ts]=tdm_init(stcfg,1,4,2,ts); % initialize model
pa.ihceq=4;     % 4=Neely synapse model
pa.hbnl=1;      % 1=nonlinear OHC transduction
            sav=tdm_step(pa,sav,cur,ts);  % step through time
wnr=sav.wnr;
end % return

% execute ecochg stimulus condition
function [ecochg,pr]=ecochg_condition(pr,inv)
ts.done = 0;
lv=pr.lv;      % tone level (dB SPL)
fr=pr.fr;      % tone frequency (kHz)
t1=pr.td;      % tone duration (msec)
t2=pr.ts;      % tone start    (msec)
tr=pr.tr;      % ramp duration (msec)
sd=pr.sd;      % ramp duration (msec)
dsp1=4;        % display #1 type
dsp2=4;        % display #2 type
if (isfield(pr,'dsp1')), dsp1=pr.dsp1;   end
if (isfield(pr,'dsp2')), dsp2=pr.dsp2;   end
if (inv), t1=-t1; end % invert stimulus
stcfg=[fr;lv;t1;t2;tr];
[pa,sav,cur,ts]=tdm_init(stcfg,1,dsp1,dsp2,ts);  % initialize model
% pin tdm25 mechanical drive for ecochg only (par_CEL25 ME tweaks + fluid coupling)
pa.rma=pa.rma/2; pa.rst=pa.rst/2;                % par_CEL25: halve malleus/stapes damping
pa.kma=pa.kma*4; pa.kst=pa.kst*4;                % par_CEL25: quadruple malleus/stapes stiffness
pa.red=pa.red*5;                                 % par_CEL25: x5 eardrum damping
pa.aflom_fac=2;                                  % tdm25 fluid coupling: cp.ac/(2*rho*dx)
pa=parm_rep(pa,pr);                              % parameter replacement
if (isfield(pr,'tsdsp')), ts.tsdsp=pr.tsdsp; end % time-step display
[sav,cur,ts,cp] = tdm_step(pa,sav,cur,ts);       % step through time
                  tdm_show(pa,sav,cur,ts,cp);    % plot saved variables
nt=length(sav.wnr);
n=nt-mod(nt,2);
tt=n*pa.dt*pa.ntsw*1000;
dt=ts.dt;
nhc=round(t2/dt/10); % why /10 ??? (matches tdm25 ecochg neural timing)
nnr=round(sd/dt)+nhc;
xl=pa.xl*10;
nx=pa.n;
xx=linspace(0,xl,nx)';
hcc=sav.hc(1:n,:)*pr.hcc; % hair-cell current
nrc=sav.nr(1:n,:)*pr.nrc; % neural-rate current
hcc(isnan(hcc))=0;
nrc(isnan(nrc))=0; % sanitize apical-boundary NaN (matches hcc); lets ihceq=0 run
hcc=[zeros(nhc,nx);hcc(1:(n-nhc),:)];
nrc=[zeros(nnr,nx);nrc(1:(n-nnr),:)];
ne=length(pr.xe);
zvc=zeros(nx,ne);
ecochg=zeros(n,ne);
for k=1:ne
    xvc=(xx-pr.xe(k))/pr.wvc;
    zvck=exp(-abs(xvc));
    ecochgk=hcc*zvck-sum(nrc,2);
    ecochgk=filter([1 -1],[1 -pr.hca],ecochgk);
    ecochg(:,k)=ecochgk;
    zvc(:,k)=zvck;
end
pr.zvc=zvc;
pr.tt=tt;
pr.xl=xl;
end % return

function p=f2p(f,a,b,c,xl)
p = (1 - (log10((f / a) + c) / b) / 100) * xl;
end

function [d1,d2,pe,ps,vr,vs,ve,nr]=fetch_sav(i,sav)
d1 = fft(sav.d1(:,i)); d2 = fft(sav.d2(:,i)); nr = fft(sav.nr(:,sav.isv(i))); vs = fft(sav.vst); ps = fft(sav.pst); pe = fft(sav.ped); ve = fft(sav.ved); vr = fft(sav.vep);
nf=length(sav.f); ii=1:nf; d1 = d1(ii); d2 = d2(ii); nr = nr(ii); vs = vs(ii); ps = ps(ii); pe = pe(ii); ve = ve(ii); vr = vr(ii);
end

function [bf,qe,mi]=find_bf(f,thr,fmin)
% APICAL ARTIFACT GUARD (2026-07-27).
% The original took a GLOBAL max over the whole frequency vector with no band
% limit and no validity test. vh carries an explicit *s = 2*pi*i*f factor, so
% the sensitivity curve is f-weighted; at apical places the response above CF is
% negligible and that weighting makes the curve rise monotonically to the TOP
% BIN, so the detector returned the band edge. Observed: the two most apical of
% the seven default places BOTH reported BF 24.98 kHz -- the last bin, identical
% to the digit -- at levels of -39 dB, against true BFs below 0.78 kHz.
%
% That single failure is the likely source of the fold==range signature that
% appeared on EVERY model this session, m=3b included: a map reading
% 24.98, 24.98, 0.78, ... has one enormous downward step spanning its full span.
%
% FIX: reject a maximum sitting on either edge of the admissible range. An edge
% maximum means the curve never turned over, so no BF was found, and NaN is the
% correct answer -- a plausible-looking wrong number is not. Callers already
% filter on isfinite. No upper band limit is imposed, so legitimate basal BFs
% near the top of the analysis range are unaffected; only a peak AT the final
% bin is rejected.
if (nargin<3 || isempty(fmin)), fmin=0.1; end
ok = f(:)>=fmin & isfinite(thr(:));
if (nnz(ok)<3), bf=NaN; qe=NaN; mi=NaN; return; end
t = thr(:); t(~ok) = -inf;
[mm,mi] = max(t);
io = find(ok);
if (mi<=io(1) || mi>=io(end))
    bf=NaN; qe=NaN; mi=NaN; return   % never turned over -> detector failed here
end
bf=f(mi); qe=bf/trapz(f,10.^((thr(:)-mm)/10));
end

function y=dbmag(x)
y=20*log10(abs(max(x,eps)));
end

function x=short_chirp(nn,n1,n2,nt)
nf=nn/2+1; mg=1; ph=0; X=zeros(nf,1);
for k=1:nf, X(k)=mg*exp(-1i*ph); dp=(2*pi)*(n1+(k/nf)*(n2-n1))/nn; ph=ph+dp; end
X=[X;zeros((nt-nn)/2,1)]; x=ffs(X); x=x-x(1); x=x*1.3183/rms(x);
end

function H=ffa(h)
H=fft(real(h)); n=length(H); m=1+n/2; H(1,:)=real(H(1,:)); H(m,:)=real(H(m,:)); H((m+1):n,:)=[];
end

function h=ffs(H)
m=length(H); n=2*(m-1); H(1,:)=real(H(1,:)); H(m,:)=real(H(m,:)); H((m+1):n,:)=conj(H((m-1):-1:2,:)); h=real(ifft(H));
end

function a=rms(x)
a=sqrt(mean(x.^2));
end

function gd=delay(R,f)
[n,m] = size(R); ph = unwrap(angle(R))/(2*pi); gd = zeros([n,m]); for k=1:m, gd(:,k) = -cdif(ph(:,k))./cdif(f(:)); end
end

function dx=cdif(x)
n=length(x); dx=zeros(size(x)); dx(1)=x(2)-x(1); dx(2:(n-1))=(x(3:n)-x(1:(n-2)))/2; dx(n)=x(n)-x(n-1);
end

function reset_color_index
if (~isoctave), ax=gca(); set(ax,'ColorOrderIndex',1); end
end

function o = isoctave
o = exist('OCTAVE_VERSION', 'builtin');
end

%==============================================================

%% modpar26, par_CEL16, modpar26c3 moved to shared files modpar26.m / modpar26c3.m

function S=opt_protocol(pr)
    if isnumeric(pr) && isscalar(pr)
        nch = max(1, pr);
    else
        nch = 1;
    end
    pa = modpar26(nch);
    fprintf('Starting optimization for nch=%d\n', nch);
    
    pv0 = [pa.ihcdr, pa.ihcrr, pa.ihcex];
    op = optimset('Display','iter','MaxIter',15);
    
    pv_opt = fminsearch(@(pv) opt_cost(pv, nch, pa), pv0, op);
    
    fprintf('Optimized parameters: ihcdr=%g, ihcrr=%g, ihcex=%g\n', pv_opt(1), pv_opt(2), pv_opt(3));
    S.pv_opt = pv_opt;
end

function err = opt_cost(pv, nch, pa)
    pa.ihcdr = pv(1);
    pa.ihcrr = pv(2);
    pa.ihcex = pv(3);
    
    % Evaluate tbabr latencies
    f=[1.0; 4.0]; slv=[40 80];
    lat=zeros(2,2); abr=zeros(2,2);
    b=12.9;c=5;d=0.413; 
    for j=1:2
        for k=1:2
            fr=f(j); lv=slv(k);
            ts.done = 0;
            td=4/sqrt(fr);
            [pa_t,sav,cur,ts]=tdm_init([fr;lv;td],nch,4,2,ts,pa);
            pa_t.ihceq=4; pa_t.hbnl=1;
            sav=tdm_step(pa_t,sav,cur,ts);
            ped=sav.ped; ped=ped - mean(ped); pmn=min(ped); pmx=max(ped);
            [~,ix]=max(sav.wnr); tpk=(ix-1)*ts.dt*ts.nw;
            lat(j,k)=tpk;
            i=lv/100; abr(j,k)=b*c^(-i)*fr^(-d);
        end
    end
    mse_tbabr = mean((lat(:) - abr(:)).^2);
    
    % Evaluate fwdmsk (amount of masking)
    amt0=[19.19 13.28; 26.46 18.31];
    mlv=[60 80]'; dly=[20 40]'; tlv=20; dt=0.02;
    prb=zeros(2,2); maxdev=6;
    % unmasked
    ts.done = 0; stcfg=[1 1;-90 tlv;100 20;0 100+dly(1)];
    [pa_t,sav,cur,ts]=tdm_init(stcfg,nch,4,2,ts,pa);
    pa_t.ihceq=4; pa_t.hbnl=1; sav=tdm_step(pa_t,sav,cur,ts);
    wnr=sav.wnr; nt=length(wnr); iprb=round((100+dly(1))/dt):nt; wnrref=max(wnr(iprb));
    for j=1:2
        for k=1:2
            plvs=tlv+amt0(j,k)+[-1 0 1]*maxdev;
            wnrmax=zeros(size(plvs));
            for i=1:length(plvs)
                ts.done = 0; stcfg=[1 1;mlv(j) plvs(i);100 20;0 100+dly(k)];
                [pa_t,sav,cur,ts]=tdm_init(stcfg,nch,4,2,ts,pa);
                pa_t.ihceq=4; pa_t.hbnl=1; sav=tdm_step(pa_t,sav,cur,ts);
                wnr=sav.wnr; wnrmax(i)=max(wnr(round((100+dly(k))/dt):length(wnr)));
            end
            if (min(wnrmax)>=wnrref), prblev=min(plvs);
            elseif (max(wnrmax)<=wnrref), prblev=max(plvs);
            else
                for i=2:length(plvs)
                    if ((wnrmax(i-1)<=wnrref)&&(wnrmax(i)>wnrref))
                        a=(wnrmax(i)-wnrref)/(wnrmax(i)-wnrmax(i-1));
                        prblev=a*plvs(i-1)+(1-a)*plvs(i); break;
                    end
                end
            end
            prb(j,k)=prblev;
        end
    end
    amt=prb-tlv;
    mse_fwdmsk = mean((amt(:) - amt0(:)).^2);
    
    err = mse_tbabr + mse_fwdmsk;
    fprintf('  MSE: ABR=%.3f, FwdMsk=%.3f, Total=%.3f\n', mse_tbabr, mse_fwdmsk, err);
end

function dof = resolve_dof(pa)
% RESOLVE_DOF  Partition DOF count. Stage 3 of the macro/micro separation.
%
% The capstone design note requires chamber count and DOF count to be
% "parameters, not forks". This was previously derived from pa.m at TWO sites
% (tdm_init and coupeig) that had to be kept in sync by hand, with a comment
% saying so. One function now serves both.
%
% DERIVED FLOOR: what the chamber configuration structurally requires.
%   2   d1 (BM) + d2 (TM-RL shear)
%   3   m>=4 (CL is a fluid chamber, d3 couples it) or pa.d3int (d3 internal)
%   4   + the resonant vent, which carries its own state. A pure-inertia vent
%       (clvk=clvr=0) stays algebraically eliminated in the a2 stamp and needs
%       no state, so dof stays 3 and prior results are untouched.
%
% EXPLICIT REQUEST may RAISE the count, never lower it. pa.dof=3 on a 1-chamber
% macromechanics is the capstone workhorse ("run 1-chamber with 3-DOF") and is
% honoured; a request BELOW the floor is ignored, because m>=4 structurally
% needs d3 to couple the CL chamber and honouring a stale lower value would
% silently produce a different model. This matters because tdm_init WRITES
% pa.dof back (pa.dof=dof), so every saved fit carries one as a side effect
% rather than as an intent -- max() makes a stale value harmless.
dof = 2;
if (pa.m>=4 || (isfield(pa,'d3int') && pa.d3int)), dof = 3; end
if (dof==3 && isfield(pa,'clvent') && pa.clvent > 0 && ...
    ((isfield(pa,'clvk') && any(pa.clvk~=0)) || (isfield(pa,'clvr') && any(pa.clvr~=0)) || ...
     (isfield(pa,'clvoct') && isfinite(pa.clvoct))))
    dof = 4;
end
if (isfield(pa,'dof') && ~isempty(pa.dof) && isfinite(pa.dof))
    dof = max(dof, round(pa.dof));
end
end
