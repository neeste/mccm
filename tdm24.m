% tdm24 - multichamber time-domain cochlear model.
%   S = tdm24(STCFG,NCH,DSP1,DSP2) computes a cochlear-model response
%   with NCH fluid chambers and stimulus described by STCFG.
%   When NCH=0, only the middle-ear stage of the model is computed.
%   STCCFG specifies the type of stimulus, for example:
%       0        - selects a click,
%       F        - selects a tone with frequency = F kHz,
%     'file.wav' - selects a waveform stimulus.
%   Otherwise, STCFG may select an experimental prtocol:
%     'tbabr'   - whole-nerve-response (WNR) latency to tonebursts,
%     'fwdmsk'  - tone-on-tone forward masking of WNR.
%     'ldngrw'  - loundess growth of single tone,
%     'dpoae'   - distortion-product otoacoustic emission (NYI).
%     'ecog'    - ECoG simulation
%   DSP1 = upper subplot [1=Vr 2=Pe 3=Ps 4=WNR]
%   DSP2 = lower subplot [1=HB 2=NR 3=LD]

function S = tdm24(stcfg,nch,dsp1,dsp2)
if (nargin<1), stcfg = 0; end
if (nargin<2), nch = 0;   end
if (nargin<3), dsp1 = 0;  end
if (nargin<4), dsp2 = 0;  end
ts = check_protocol(stcfg,nch,dsp1,dsp2); % protocol ?
if (ts.done), if (nargout), S = ts; end; return; end
if (nargin<3), dsp1 = 1;  end
if (nargin<4), dsp2 = 2;  end
% compute time-domain model response ------------
[pa,sav,cur,ts] = tdm_init(stcfg,nch,dsp1,dsp2,ts); % initialize model
[sav,cur,ts,cp] = tdm_step(pa,sav,cur,ts);          % step through time
                  tdm_show(pa,sav,cur,ts,cp);       % plot saved variables
% -----------------------------------------------
if (nargout), S = sav; end    % return structure of saved variables
end % return

%==========================================================
% initialize stimulus, model variables, and state
function [pa,sav,cur,ts]=tdm_init(stcfg,nch,dsp1,dsp2,ts)
% specify mechanical parameters
pa=modpar24(nch);
if (nch<1)       % cochlea impedance only
    pa.ihcv = 0; %
end
if (pa.nmev==1) % simplified middle-ear
    if (pa.stim==3)
        pa.stim=2;
        fprintf('moved earphone stimulus to eardrum\n')
    end
end
if(~isfield(pa,'ihcxx')), pa.ihcxx=[]; end
if(~isfield(ts,'ohcxx')), ts.ohcxx=[]; end
if (length(ts.ohcxx)>1)
    if (ts.ohcxx(2)>ts.ohcxx(1))
        nn=1+round((pa.n-1)*(ts.ohcxx/(10*pa.xl)));
        ts.ohcxx=false(pa.n,1);
        ts.ohcxx(nn(1):nn(2))=true;
    else
        ts.ohcxx=[];
    end
end
% specify timing parameters
pa.nstp = 20480;            % number of time steps per cycle
pa.ncyc = 1;                % number of cycles
pa.ntsw = 10;               % number of time steps between saves
pa.ntsf = 2048;             % number of time steps saved
pa.dt = 2e-6;               % time step (s)
% set numerical integration parameter
pa.nimp = 1;                % number of "improvements"
% initialize state
dof = 2;                    % degress of freedom
pa.ncp = pa.n * dof;        % number of BM variables
if (pa.nmev<1), pa.nmev = 1; end % need at least 1 ME variable
nsv = pa.ncp + pa.nmev;      % number of state variables
cur.d = zeros(nsv,1);        % displacement
cur.v = zeros(nsv,1);        % velocity
cur.hb = zeros(pa.n,1);      % HB displacement
cur.nr = zeros(pa.n,1);      % neural rate
cur.ld = zeros(pa.ldne,1);  % partial loudness
cur.stp = 0;
cur.cyc = 0;
cur.wnr = 0;
cur.vst = 0;
cur.pst = 0;
cur.ped = 0;
cur.ved = 0;
cur.vep = 0;
cur.stm = 0;
cur.nrme = zeros(pa.nmev,1);
% initialize saved variables
nt = pa.ntsf;
nsv = length(pa.isv);
sav.vep = zeros(nt,1);
sav.ped = zeros(nt,1);
sav.ved = zeros(nt,1);
sav.pst = zeros(nt,1);
sav.vst = zeros(nt,1);
sav.wnr = zeros(nt,1);
sav.d1 = zeros(nt,nsv);
sav.d2 = zeros(nt,nsv);
sav.hb = zeros(nt,pa.n);
sav.nr = zeros(nt,pa.n);
sav.ld = zeros(nt,pa.ldne);
sav.nch = nch; % number of fluid chambers
% initialize stimulus & graphics
ts=ts_init(pa,sav,cur,ts,1,2,stcfg,dsp1,dsp2); % figure numbers, total time
if (ts.tt == 0), return; end            % no stimulus
ts.smx = 0;
% initialize frequency
nw = ts.nw;
nt = ts.nt / nw;
dt = ts.dt * nw / 1000;
sav.f = (0:((nt/2) - 1))' / (1000 * dt * nt);
sav.smx = ts.smx;
end % return

% step through time
function [sav,cur,ts,cp]=tdm_step(pa,sav,cur,ts)
cp=[];
if (pa.m<0), ts=plot_spec(pa,sav,ts); end % spectrogram
if (ts.done), return; end                 % done
% prepare BM & fluid equations -----------------------
      me=midear(pa);     % middle-ear % coupler
 [cp,me]=cochlea(pa,me); % cochlear partition & fluid
[cur,nf]=neural(pa,cur); % sensory cell, synaptic, & neural
if (~isempty(ts.ohcxx)), cp.gm(ts.ohcxx)=0; end
% step through time ----------------------------------
nxt = applst(pa,ts,cur);                % apply stimulus to initial state
sav = savest(pa,sav,cur,ts);            % save initial state
while (cur.cyc < pa.ncyc)               % time-step loop ---------------
    nxt = inctim(cur,nxt,ts);           % increment state to first step
  	[cur,ac] = accel(pa,cp,me,cur,cur); % get current acceleration
    cur = ihc_step(pa,cur,nf);          % update HB state
    sav = savest(pa,sav,cur,ts);        % save current state
  	nxt = update(pa,cur,nxt,ac,ac);     % integrate to next step
    nxt = applst(pa,ts,nxt);            % apply stimulus to next state
  	for imp=1:pa.nimp                   % improve integration  --
		[cur,an] = accel(pa,cp,me,cur,nxt); % get next acceleration
     	nxt = update(pa,cur,nxt,ac,an); % integrate state to next step
    end                                 % end improvement loop --
	cur=nxt;                            % make next state current
end                                     % end time-step loop -----------
end % return

%==========================================================

function me=midear(pa)
% set all middle-ear variables to zero
nme=pa.nmev; % number of middle-ear components
me.pedm=zeros(1,nme); me.pedr=zeros(1,nme); me.pedk=zeros(1,nme);
me.smes=zeros(nme,1); me.qstr=zeros(1,nme); me.qstk=zeros(1,nme);
me.smem=zeros(nme,nme); me.smer=zeros(nme,nme); me.smek=zeros(nme,nme);
if (pa.nmev==1) % simplified middle-ear
    me.smem = 1;
    me.smer = (pa.rme / pa.mme);
    me.smek = (pa.kme / pa.mme);
    me.qstr = pa.rme;
    me.qstk = pa.kme;
    me.smes = -1;
    me.ied = 1; me.ist=1;
    me.idx = [1 1]; % ME indices for eardrum & stapes
    me.vst = [0 0 0 1];
    return
end
% middle-ear-component indices
iep=1; % earphone diaphragm
icp=2; % eardrum (non-propagating motion)
ima=3; % malleus
ist=4; %s stapes
me.idx=[icp ist]; % ME indices for eardrum & stapes
% middle-ear useful constructs
cr = pa.cep;
ad = pa.adi;
ac = pa.acp;
au = pa.aed;
am = pa.ama;
as = pa.ast;
gm = pa.gme;
md = pa.mdi;
mc = pa.mcp;
mu = pa.med;
mm = pa.mma;
ms = pa.mst;
ru = pa.red;
rm = pa.rma;
ri = pa.rim;
rs = pa.rst;
rc = pa.rcp;
ra = pa.rfz;
rd = pa.rdi + pa.rep;
ku = pa.ked;
km = pa.kma;
ki = pa.kim;
ks = pa.kst;
ka = pa.kcp;
kd = pa.kdi;
% set middle-ear mass
me.smem(iep,iep) = md / ad^2;
me.smem(icp,icp) = mc / ac^2 + mu  / au^2;
me.smem(icp,ima) = -mu  / au^2;
me.smem(ima,icp) = -mu  / au^2;
me.smem(ima,ima) = mm / am^2 + mu  / au^2;
me.smem(ist,ist) = ms / as^2;
% set middle-ear damping
me.smer(iep,iep) = (rd + ra) / ad^2;
me.smer(iep,icp) = -ra / (ad * ac);
me.smer(icp,iep) = -ra / (ad * ac);
me.smer(icp,icp) = (rc + ra) / ac^2 + ru / au^2;
me.smer(icp,ima) = -ru / au^2;
me.smer(ima,icp) = -ru / au^2;
me.smer(ima,ima) = (rm +  gm^2 * ri) / am^2 + ru / au^2;
me.smer(ima,ist) = -gm * ri / (am * as);
me.smer(ist,ima) = -gm * ri / (am * as);
me.smer(ist,ist) = (rs + ri) / as^2;
% set middle-ear stiffness
me.smek(iep,iep) = (kd + ka) / ad^2;
me.smek(iep,icp) = -ka / (ad * ac);
me.smek(icp,iep) = -ka / (ad * ac);
me.smek(icp,icp) = ka / ac^2 + ku / au^2;
me.smek(icp,ima) = -ku / au^2;
me.smek(ima,icp) = -ku / au^2;
me.smek(ima,ima) = (km +  gm^2 * ki) / am^2 + ku / au^2;
me.smek(ima,ist) = -(gm * ki) / (am * as);
me.smek(ist,ima) = -(gm * ki) / (am * as);
me.smek(ist,ist) = (ks + ki) / as^2;
% invert mass matrix
me.smem=inv(me.smem);
% apply stimulus
if (pa.stim==3)     % --> earphone
    me.smes(iep) = -(cr / ad);
elseif (pa.stim==2) % -->  eardrum
    me.smes(icp) = -1;
    me.smes(ima) = -1;
end
% pressure recovery
me.pedm(icp) = mu / au^2;
me.pedr(icp) = ru / au^2;
me.pedk(icp) = ku / au^2;
me.pedm(ima) = -mu / au^2;
me.pedr(ima) = -ru / au^2;
me.pedk(ima) = -ku / au^2;
% stapes boundary condition
qstrst = (rs + ri) / ms;
qstkst = (ks + ki) / ms;
qstrma = (gm * as * ri) / (am * ms);
qstkma = (gm * as * ki) / (am * ms);
if (pa.m)
    mst = ms / as^2;
    me.qstr(ist) = -mst * qstrst;
    me.qstk(ist) = -mst * qstkst;
    me.qstr(ima) =  mst * qstrma;
    me.qstk(ima) =  mst * qstkma;
else
    asco = as + pa.mco * as^2 / ms;
    me.qstr(ist) = (pa.rco - pa.mco * qstrst) / asco;
    me.qstk(ist) = (       - pa.mco * qstkst) / asco;
    me.qstr(ima) = (         pa.mco * qstrma) / asco;
    me.qstk(ima) = (         pa.mco * qstkma) / asco;
    me.mem = 0;
end
me.mst = (ms / as);
me.ved = [0 10 0 0]; % <-- why 10 ??
if (pa.m==1)
    me.vst = [0 0 0 2]; % why 2 ??
else
    me.vst = [0 0 0 1];
end
end % return

function [cp,me]=cochlea(pa,me)
m = pa.m;
n = pa.n;
cp = imped(pa);
if (m<3) % specify HB displacement
    cp.hb=[cp.gh -ones(n,1)];
elseif (m==3)
    cp.hb=[zeros(n,1) ones(n,1)];
end
dx = pa.xl / (n - 1);                    % radial slice thickness
if (m < 1) % no fluid chambers ??
    cp.abmom=0;
    cp.alfx = me.mst / (2 * pa.rho * dx);
    return;
end
% setup m-chamber fluid equations
cp.abmom = (cp.bw * dx) ./ cp.m1;          % BM area over mass
cp.alfx = (pa.ast / pa.mst) / cp.abmom(1); % stapes area over mass
aflom = cp.ac / (2 * pa.rho * dx);         % fluid area over mass
kk = 2:(n-1);
mm = m*m;
a1=zeros(n,mm);a2=zeros(n,mm);a3=zeros(n,mm);
if (m<3)
    a3(1,1) = -aflom(1) ./ cp.abmom(1);      % basal BC: ME impedance
    a2(1,1) =  1 + cp.alfx - a3(1);          % "
    a1(kk,1) = -aflom(kk-1) ./ cp.abmom(kk); % left of diagonal
    a3(kk,1) = -aflom(kk) ./ cp.abmom(kk);   % right of diagonal
    a2(kk,1) = 1 - a1(kk) - a3(kk);          % diagonal element
    a1(n,1) = -1;                            % apical BC: p(n)-p(n-1)=q(n)
    a2(n,1) = 1;                             % "
else
    error('m>=3 not yet implemented');
end
if (m>1), a2(:,(1:m)+(mm-m)) = 1; end         % enforce fluid volume
% expand matrix A to multiple chambers
cp.aa=xpnd_a(a1,a2,a3,m,n);
end % return

function [cur,nf]=neural(pa,cur)
n = pa.n;
% initialize OHC
cur.vohc = zeros(n,1);    % OHC voltage
cur.dvohc = zeros(n,1);   % OHC voltage derivative
% initialize IHC
cur.vi = zeros(n,1);      % IHC receptor voltage
cur.cc = zeros(n,1);      % IHC receptor charge
cur.ss  = ones(n,1);      % IHC synapse reservoir
cur.nr  = zeros(n,1);     % IHC synapse neural rate
% IHC & synapse
nf.nihc = 10;      % IHC model time interval
nf.sf = pa.ihcsf;  % IHC membrane-potential scale factor
% model time step adjustment
dt = pa.dt * nf.nihc; % IHC model time step
nf.kk = exp(-dt / pa.ihctc);
nf.rr = pa.ihcrr * (dt / 2e-6)^(1/pa.ihcex);
nf.dr = pa.ihcdr * (dt / 2e-6)^(1/pa.ihcex);
% initialize neural ensemble
ehw=floor(pa.n/(pa.xl*2/pa.ew)); % ensemble half-width
i1k = pa.isv(3);               % 1 kHz place
nf.ens = (i1k - ehw) : (i1k + ehw);
end % return

%-------------------------------------------

function aa=xpnd_a(a1,a2,a3,m,n)
nm = n*m;
nd = 1+2*m;
ad = zeros(nm,nd);
dd = zeros(nm,nd);
% transform m-blocks to nd-diagonals
for k=1:n
    for j=1:m
        jj = (1:m)+(j-1)*m;
        ii = (1:nd)+(j-1);
        kk = j+(k-1)*m;
        aaa = [a1(k,jj) a2(k,jj) a3(k,jj)];
        ad(kk,:) = aaa(ii);
    end
end
% skew diagonals for spdiags
for k=1:nm
    for j=1:m
        dd(k,1+m)                   = ad(k,1+m);
        if (k>j),      dd(k-j,1+m-j) = ad(k,1+m-j); end
        if (k<=(nm-j)), dd(k+j,1+m+j) = ad(k,1+m+j); end
    end
end
% format as sparse multi-diagonal matrix
aa = spdiags(dd,-m:m,nm,nm);
end % return

function qx=xpnd_q(qst,ss,m,n)
s1 = ss(:,1);
qst = [qst;zeros(n-1,1)];
qx = zeros(n*m,1);
if (m==1)
    qx = s1 + qst;
elseif (m==2)
    j2 = (1:n) * 2;
    j1 = j2 - 1;
    qx(j1) = s1 + qst;
elseif (m==3)
    j3 = (1:n) * 3;
    j2 = j3 - 1;
    j1 = j2 - 1;
    qx(j1) = s1 + qst;
end
end % return

%-------------------------------------------
%
% plot saved variables
%
function tdm_show(pa,sav,cur,ts,cp)
report_level(pa,sav,ts); % compute SPL at eardrum
plot_ts(pa,sav,cur,ts);  % time-step plot (summary)
plot_sv(pa,sav,ts,cp);   % saved-variables plot
end % return

%==========================================================

% distribute cochlear-partition properties over x
function cp=imped(pa)
n = pa.n;
dx = pa.xl / (n - 1);
x = (0:(n - 1))' * dx;
x = x.*(1+(pa.xtap*x).^pa.xtex);
q = x .^ 2;
cp.k1 = pa.k1o * exp(pa.k1e * x + pa.k1q * q);
cp.r1 = pa.r1o * exp(pa.r1e * x + pa.r1q * q);
cp.m1 = pa.m1o * exp(pa.m1e * x + pa.m1q * q);
cp.k2 = pa.k2o * exp(pa.k2e * x + pa.k2q * q);
cp.r2 = pa.r2o * exp(pa.r2e * x + pa.r2q * q);
cp.m2 = pa.m2o * exp(pa.m2e * x + pa.m2q * q);
cp.k3 = pa.k3o * exp(pa.k3e * x + pa.k3q * q);
cp.r3 = pa.r3o * exp(pa.r3e * x + pa.r3q * q);
cp.k4 = pa.k4o * exp(pa.k4e * x + pa.k4q * q);
cp.r4 = pa.r4o * exp(pa.r4e * x + pa.r4q * q);
cp.gh = pa.gpo * exp(pa.gpe * x + pa.gpq * q);
cp.ac = pa.aco * exp(pa.ace * x + pa.acq * q);
cp.bw = pa.bwo * exp(pa.bwe * x + pa.bwq * q);
cp.gm = pa.gam * ones(size(x));
% apical BC
ihe=n;
cp.k1(ihe) = pa.khe;
cp.r1(ihe) = pa.rhe;
cp.m1(ihe) = pa.mhe;
cp.k2(ihe) = 0;
cp.r2(ihe) = 0;
cp.m2(ihe) = 0;
cp.k3(ihe) = 0;
cp.r3(ihe) = 0;
cp.k4(ihe) = 0;
cp.r4(ihe) = 0;
end % return

% increment time
function nxt=inctim(cur,nxt,ts)
nxt.stp = cur.stp + 1;
nxt.cyc = cur.cyc;
if (nxt.stp >= ts.nt )
   nxt.stp = nxt.stp - ts.nt;
   nxt.cyc = nxt.cyc + 1;
end
end % return
%-------------------------------------------

% apply stimulus to given state
function st=applst(pa,ts,st)
if (isempty(ts.stwav)), return; end
if (st.stp>=ts.nt)
    stm = 0;
else
    stm = ts.stwav(st.stp+1);
end
switch (pa.stim)
    case 1, st.pst=stm;
    case 2, st.ped=stm;
    case 3, st.vep=stm;
end
st.stm=stm;
st.qme(pa.stim)=stm;
end % return

%-------------------------------------------
% determine acceleration given state
function [cur,a]=accel(pa,cp,me,cur,st)
a = zeros(size(st.d));
stm = st.stm;
% ME pressure
nme = pa.nmev;
ime = pa.ncp + (1:nme);
d0 = st.d(ime);
v0 = st.v(ime);
ist = me.idx(2); % ME index for stapes
% BM pressure
[ss,ii]=force_cp(pa,cp,st);
% micromechanics
s0 = stm * me.smes - (me.smek * d0 + me.smer * v0);
prw = (pa.rrw * v0(ist) + pa.krw * d0(ist)) * pa.arw;
m = pa.m;
n = pa.n;
if (pa.stim==1) % stapes pressure stimulus
    qst = stm;
elseif (m)
    qst = -cp.alfx * (s0(ist) + prw);
else
    qst = me.qstr * v0 + me.qstk * d0;
end
if (m)
    qq=xpnd_q(qst,ss,m,n); % expand to multiple chambers
    pp = cp.aa \ qq;       % solve matrix equation
else
    pp = qst;            % cochlea input impedance only
end
% ME acceleration
s0 = stm * me.smes - (me.smek * d0 + me.smer * v0);
s0(ist) =  s0(ist) + pp(1) + prw;
a0 = me.smem * s0;
a(ime)=a0;      % return middle-ear acceleration
% cochlea acceleration
if (m)
    a = fold_p(pp,m,n,a,ss,ii,cp);
else                       % phantom cochlea
    a(ii(:,1)) = -a0(ist); % set BM_acceleration = stapes_acceleration
end
% current pressure
if (pa.stim==2) % copy eardrum pressure to receiver voltage
    cur.vep = cur.ped;
else            % update eardrum pressure
    cur.ped = me.pedm * a0 + me.pedr * v0 + me.pedk * d0;
end
cur.ved = me.ved * v0;
cur.vst = me.vst * v0;
cur.pst = pp(1);
st.a = a;
end % return
%-------------------------------------------
function a=fold_p(pp,m,n,a,ss,ii,cp)
s1=ss(:,1);
s2=ss(:,2);
i1=ii(:,1);
i2=ii(:,2);
if (m==1)
    a(i1) = (s1 - pp) ./ cp.m1;
    a(i2) = s2 ./ cp.m2;
elseif (m==2)
    j2=(1:n)*2;
    j1=j2-1;
    s1 = s1 - (pp(j1) - pp(j2)) / 2; % (ST-SV) pressure difference
    a(i1) = s1 ./ cp.m1;
    a(i2) = s2 ./ cp.m2;
elseif (m==3)
    j3=(1:n)*3;
    j2=j3-1;
    j1=j2-1;
    s1 = s1-(pp(j1) - pp(j3)); % (ST-SV) pressure difference
    s2 = s2-(pp(j2) - pp(j3)); % (SS-SV) pressure difference
    a(i1) = s1 ./ cp.m1;
    a(i2) = s1 ./ cp.m1 + s2 ./ cp.m2;
end
a(i1(1)) = 0;
a(i1(n)) = 0;
end % return
%-------------------------------------------
% integrate state to next time step
function nxt=update(pa,cur,nxt,cur_a,nxt_a)
dt2 = pa.dt / 2;
nxt.v = cur.v + (cur_a + nxt_a) * dt2;
nxt.d = cur.d + (cur.v + nxt.v) * dt2;
% IHC
nxt.vi = cur.vi;
nxt.ss = cur.ss;
nxt.nr = cur.nr;
nxt.ld = cur.ld;
end % return

% update HB state
function cur=ihc_step(pa,cur,nf)
n = pa.n;
ii=1:n;
if (pa.ihcv==0) % HB displacement or velocity ?
    d1 = cur.d(ii);
    d2 = cur.d(ii+n);
    %cur.hb = cp.hb(:,1) .* d1 + cp.hb(:,2) .* d2; % HB displacement
    cur.hb = d1 - d2; % HB displacement
else
    v1 = cur.v(ii);
    v2 = cur.v(ii+n);
    %cur.hb = cp.hb(:,1) .* v1 + cp.hb(:,2) .* v2; % HB velocity
    cur.hb = v1 - v2; % HB velocity
end
if (pa.ihceq<1), return; end
if ((nf.nihc>1) && mod(cur.stp,nf.nihc)), return; end
n = pa.n;
% synapse
ex = pa.ihcex;  % IHC synapse exponent
% model time step adjustment
dt = pa.dt * nf.nihc; % IHC model time step
sf = nf.sf;  % IHC membrane-potential scale factor
kk = nf.kk;  % first-order LPF
rr = nf.rr;  % reservoir recovery rate
dr = nf.dr;  % reservoir depletion rate
% loop across all IHCs
wncc = 0;
for i=1:n
    % IHC transduction
    s = cur.hb(i);                     % HB displacement
    s = max(0,s);                      % rectify
    s = (1 - kk) * s + kk * cur.vi(i); % LPF
    cur.vi(i) = s;                     % save IHC membrane potential
    % IHC synapse
    vi = s * sf;                       % scale IHC membrane potential
    ss = cur.ss(i);                    % reservoir contents
    nr = cur.nr(i);                    % neural rate
    ss = ss + (rr * (1 - ss) - dr * nr) * dt;
    ss = max(0, min(ss, 1));
    nr = ss^ex * vi;
    cur.ss(i) = ss;
    cur.nr(i) = nr;
    wncc = wncc + cur.nr(i);           % whole-nerve cumulative count
end
if (~isempty(pa.ihcxx)), cur.nr(pa.ihcxx)=0; end
cur.wnr = sum(cur.nr) * pa.nrgn;       % scale whole-nerve sum of rates
ew = pa.ldew/(pa.xl*10);
ne = pa.ldne;
sc = pa.ldsc;
cur.ld = soft_max(loudness(cur.nr,ew,ne,sc).^(1/3),10);
end % return

function x=soft_max(x,m)
x=x./(m+sum(x));
end % return

% cochlear-partition force
function [ss,ii]=force_cp(pa,cp,st)
n = pa.n;
i1 = 1:n;
i2 = (n+1):(2*n);
ii=[i1(:) i2(:)];
if (pa.m<1), ss=0; return; end
d1 = st.d(i1);
v1 = st.v(i1);
d2 = st.d(i2);
v2 = st.v(i2);
gam = cp.gm;
if (pa.hbnl)
    d3 = cp.hb(:,1) .* d1 + cp.hb(:,2) .* d2;
    if (pa.mmeq == 1)          % MOH90
        hbt = max(abs(d3) / pa.hbmx,1);
        gam = cp.gm ./ (1 + pa.hbsc * log(hbt));
    elseif (pa.mmeq == 9)      % CEL16 & MOH17
        dbt = abs(d3) / pa.hbmx;
        if (dbt > 1)
            gam = cp.gm / (1 + pa.hbsc * log(dbt));
        end
    end
end
if (pa.m<3)
    d3 = d1 - d2;
    v3 = v1 - v2;
    s1tmp = cp.k1 .* d1 + cp.r1 .* v1;
    s2tmp = cp.k2 .* d2 + cp.r2 .* v2;
    s3tmp = cp.k3 .* d3 + cp.r3 .* v3;
    s4tmp = cp.k4 .* d3 + cp.r4 .* v3;
    s1 = -(s1tmp + s3tmp .* cp.gh - s4tmp .* gam);
    s2 = -(s2tmp - s3tmp);
elseif (pa.m==3)
    s11 =  cp.r1           .* v1 + cp.k1                          .* d1;
    s12 = (cp.gh .* cp.r3) .* v2 + (cp.gh .* cp.k3 - gam .* cp.k4).* d2;
    s21 = cp.r2            .* v1 + cp.k2                          .* d1;
    s22 = (cp.r2 + cp.r3)  .* v2 + (cp.k2 + cp.k3)                .* d2;
    s1 = -(s11 + s12);
    s2 = -(s22 - s21);
end
s1(n) = 0;
s2(n) = 0;
ss=[s1 s2];
end % return

% save data from current state
function sav=savest(pa,sav,cur,ts)
istp = cur.stp;
n=pa.n;
nw=ts.nw;
kw=mod(istp,nw);
if (kw == 0)
    i = 1 + mod(round(istp / nw), ts.nt);
   sav.vep(i) = cur.vep;
   sav.ped(i) = cur.ped;
   sav.ved(i) = cur.ved;
   sav.pst(i) = cur.pst;
   sav.vst(i) = cur.vst;
   sav.wnr(i) = cur.wnr;
   sav.d1(i,:) = cur.d(pa.isv);
   sav.d2(i,:) = cur.d(pa.isv+n);
   sav.hb(i,:) = cur.hb;
   sav.nr(i,:) = cur.nr;
   sav.ld(i,:) = cur.ld;
   % time-step display
   ts.smx = sav.smx;
   ts = plot_ts(pa,sav,cur,ts,i);
   sav.smx = ts.smx;
end
end % return

%==========================================================

% initialize time-step + saved-state displays
function ts=ts_init(pa,sav,cur,ts,tsfig,svfig,stcfg,dsp1,dsp2)
if (~dsp1 && ~dsp2), tsfig = 0; end
ts.tsdsp = tsfig;             % time-step figure number
ts.svdsp = svfig;             % saved-variables figure number
ts.stcfg = stcfg;             % stimulus configuratopm
ts.lab = sprintf('nch=%d',pa.m);
[pa,stwav,tt,nt,smx]=stim_init(pa,stcfg); % initialize stimulus
ts.stwav = stwav; % stimulus waveform
ts.dsp1 = dsp1;   % 1=Vr 2=Pe 3=Ps 4=wnr
ts.dsp2 = dsp2;   % 1=HB 2=NR
ts = plot_ts(pa,sav,cur,ts,tt,smx); % initialize time-step plot
ts.tt = tt;                         % time length (sec)
ts.nt = nt;                         % time length (samples)
ts.dt = tt / nt;                    % time-step
ts.nw = pa.ntsw;                    % save interval
ts.ntsf = floor(ts.nt/ts.nw);       % number of time steps saved
end % return

function [pa,stwav,tt,nt,smx]=stim_init(pa,stcfg)
% stcfg has several formats to describe a variety of stimulus types.
%
% When stcfg is a character array ending with '.wav' then the stimulus
% is read from the specified waveform file. The default stimulus level
% is 60 dB SPL at the eardrum.
%
% When stcfg is a character array without '.wav' then may specify a
% an experimental protocol, such as 'tbabr' or 'fwdmsk'.
%
% When stcfg is a cell array, then the first member may specify a
% waveform and the second member its level.
%
% When stcfg=0, the stimulus is click.
%
% When stcfg is a numeric array with one row specifies a combination
% of tones at the specified frequencies (kHz). The default level is
% 40 dB SPL at the eardrum.
%
% When stcfg is a numeric array with two rows, the second row specifies
% the indiviual levels of the tone frequencies listed in the first row.
%
% When stcfg is a numeric array with three rows, the stimulus is a
% combination of Blackman-windowed tone bursts at frequencies and levels
% specified by the first two rows, and durations (msec) specified on the
% third row.
%
% When stcfg is a numeric array with three rows, the stimulus is a
% combination of Tukey-windowed tone bursts at frequencies, levels, and
% durations (msec) specified by the first three rows, and
% start times (msec) specified on the fourth row.
%
% When stcfg is a numeric array with three rows (and a negative frequency),
% then the stimulus is a specially designed for the loudness-growth
% protocol. The second row specifies the center frequency and the third
% row specifies the bandwidth. The stimulus is either five-tone (when
% bandwidth is positive) or flat-noise (when bandwidth is negative).
%
sr=1/pa.dt;
dt=2e-5; % 20 usec
nt=pa.nstp;
nn=round(nt*pa.dt/dt);
tt=nn*dt*1000;         % total time
if (iscell(stcfg))     % waveform stimulus (scaled)
    [yy,nt,tt]=waveform(stcfg,sr);
elseif (ischar(stcfg))
    [yy,nt,tt]=waveform(stcfg,sr);
elseif (stcfg==0)      % chirp stimulus
    n1=1+round(0.001/dt);
    n2=1+round(0.002/dt);
    yy=short_chirp(nn,n1,n2,nt);
elseif (size(stcfg,1)<=2) % tone stimulus
    [yy,nt,tt]=multi_tone(stcfg,tt,nt,0);
elseif (stcfg(1,1)>0)     % toneurst stimulus
    [yy,nt,tt]=multi_toneburst(stcfg,tt,nt);
else                      % band-limit stimulus
    st=stcfg(4,1); % stimulus type: 0=tone 1=flat_noise 2=lono_noise 3=AM_tone
    rc=stcfg(5,1); % random count
    [yy,nt,tt]=band_limit(stcfg,tt,nt,st,rc);
end
stwav=yy * pa.stgain;
smx=max(max(abs(stwav)),1e-9);
end % return

% waveform stimulus
function [x,nt,tt]=waveform(stcfg,sr)
spl_ref = 0.0002; % SPL reference pressure (rms Pa)
if (iscell(stcfg))
    cfn = 'pedcal.mat';
    if (~exist(cfn,'file'))
        error('needs erdrum-pressure calibration (%s)',cfn); end
    load(cfn,'f','pc_gn')
    fn=stcfg{1};
    lv=stcfg{2};
    fr=1; % calibrate 1 kHz
    if(contains(fn,'500')), fr=0.5;
    elseif(contains(fn,'1000')), fr=1;
    elseif(contains(fn,'2000')), fr=2;
    elseif(contains(fn,'4000')), fr=4;
    end
    pr=exp(interp1(log(f+eps),abs(pc_gn),log(fr)));
    gn=(spl_ref/pr/sqrt(2))*10^((lv+40)/20);
else
    fn=stcfg;
    gn=10;
end
fprintf('stimulus file: %s\n',fn);
[x,fs]=audioread(fn);
rsr=round(sr/fs);
x=interp(x,rsr);
if (contains(fn,'tb'))
    x=x*gn*2*sqrt(2)/(max(x)-min(x)); % normalize rms average
else
    x=x*gn/rms(x); % normalize rms average
end
nt=length(x);
tt=1000*nt/sr;
end % return

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
    ts=tbabr_protocol;                                  % TBABR
elseif (contains(stcfg,'fwdmsk'))
    ts=fwdmsk_protocol(pr);                             % FWDMSK
elseif (contains(stcfg,'ecog'))
    ts=ecog_protocol(pr);                               % ECOG
elseif (contains(stcfg,'ldngrw'))
    ts=ldngrw_protocol(pr,st,rc);                       % LDNGRW
elseif (contains(stcfg,'dpoae'))
    ts=dpoae_protocol(pr);                              % DPOAE
elseif (contains(stcfg,'exproto'))
    ts=user_protocol(pr);                               % EXPROTO
end
ts.done = 1;                 % protocol completed!
end % return

% TBABR protocol
function S=tbabr_protocol
fprintf('tone-burst ABR protocol...\n')
f=(0.5*2.^(0:3))';  % stimulus frequencies (kHz)
slv=[20 40 60 80];  % stimulus levels
lev=zeros(4,4);     % eardrum levels
lat=zeros(4,4);     % WNR latencies
perform_calibration % calibrate, if necessary
tic;                % start stopwatch
% collect data
for j=1:4
    for k=1:4
        fr=f(j);
        lv=slv(k);
        [mlv,tpk] = tbabr_condition(fr,lv);
        lev(j,k)=mlv;
        lat(j,k)=tpk;
        fprintf('%5.2f %3.0f %5.1f %4.1f\n',fr,lv,mlv,tpk);
    end
end
fprintf('tbabr protocol: elapsed time = %.3f (sec)\n',toc)
% analyze data
figure(1);clf
loglog(f,lat,'bo-');
xlabel('frequency (kHz)')
ylabel('latency (ms)')
axis([0.250 8 1 30])
write_data('tbabr.txt',[f lev lat]);
S.f=f; S.lev=lev; S.lat=lat;
end % return

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

% ECoG  protocol
function S=ecog_protocol(pr)
fprintf('ECoG protocol...\n')
ncnd=length(pr.lv);  % number of conditions
if (~isfield(pr,'fr'))   , pr.fr=500; end
if (~isfield(pr,'lv'))   , pr.lv=60; end
if (~isfield(pr,'td'))   , pr.td=20; end
if (~isfield(pr,'tr'))   , pr.tr=4; end
if (~isfield(pr,'ts'))   , pr.ts=0; end
if (~isfield(pr,'tt'))   , pr.tt=40; end
if (~isfield(pr,'wd'))   , pr.wd=1; end
if (~isfield(pr,'ihcxx')), pr.ihcxx=zeros(ncnd,2); end
if (~isfield(pr,'ohcxx')), pr.ohcxx=zeros(ncnd,2); end
col='rgbk';
perform_calibration % calibrate, if necessary
tic;                % start stopwatch
for k=1:ncnd
    wnr0=ecog_condition(pr,k,0);
    wnr1=ecog_condition(pr,k,1);
    if (k==1)
        nt=length(wnr0);
        wnr=zeros(nt,ncnd,2);
    end
    wnr(:,k,1)= wnr0;
    wnr(:,k,2)= wnr1;
    ihcxx=pr.ihcxx(k,:);
    ohcxx=pr.ohcxx(k,:);
end
fprintf('ecog protocol: elapsed time = %.3f (sec)\n',toc)
% plot results
lab='';
if (ihcxx(2)>ihcxx(1)), lab=[lab 'IHC ']; end
if (ohcxx(2)>ohcxx(1)), lab=[lab 'OHC ']; end
r11=wnr(:,1,1);
r12=wnr(:,1,2);
r21=wnr(:,2,1);
r22=wnr(:,2,2);
r1d=(r11-r12)/2;
r2d=(r21-r22)/2;
nt=length(r11);
t=linspace(0,pr.tt,nt)';
figure(2);clf       % clear figure 2
subplot(2,1,1)
plot(t,r11,col(1),t,r12,col(2)); hold on
plot(t,r21,col(3),t,r22,col(4)); hold on
legend('cond1,ph1','cond1,ph2','cond2,ph1','cond2,ph2','Location','northeast')
xlabel('time (ms)')
ylabel('WNR')
title(lab)
subplot(2,1,2)
plot(t,r1d,col(1),t,r2d,col(3))
legend('cond1','cond2','Location','northeast')
xlabel('time (ms)')
ylabel('ECoG')
% save wnr data
S.mlv=pr.lv;
if (pr.wd)
    wnr1=squeeze(wnr(:,:,1));
    wnr2=squeeze(wnr(:,:,2));
    wnrd=(wnr1-wnr2)/2;
    write_data('ecog.txt',[t wnr1 wnr2 wnrd]);
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
pa=modpar24(1);
fprintf('ensemble: ldew=%.1f; ldne=%d\n',pa.ldew,pa.ldne);
fprintf('  fr   lv   bw   ldph\n');
tic;                % start stopwatch
% collect data
for j=1:nbw
    bw=sbw(j);
    dn=sprintf('bw%03d',round(abs(bw),2)*100);
    savnr=exist(dn,'dir')&&(st==1);
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
            fn=sprintf('%s/ldngrw%03d.mat',dn,lv);
            nr=sav.nr;
            fr=sfr;
            save(fn,'nr','fr','bw','bwk');
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
    legend('tdm24','Fletcher','Location','SouthEast')
    axis([-5 105 0.001 100])
    subplot(1,2,2)
    plot(slv,lph,lvd,ldp,'k:')
    xlabel('level (dB SPL)')
    ylabel('loudness (phon)')
    axis([-5 105 -5 105])
    legend('tdm24','Rasetshwane','Location','NorthWest')
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

function [pld,enr]=loudness(nr,ewn,ne,c)
nx=size(nr,1);
kw=round(nx*ewn)-1;
enr=zeros(ne,1);
p=1/3; % log-loudness slope
for k=1:ne
    k1=1+round((k-1)*(nx-1)/(ne-1));
    k2=min(k1+kw,nx);
    kk=k1:k2;
    enr(k)=mean(mean(abs(nr(kk,:)).^p));
end
N=enr.^(1/p);
pld=(N.^c(1)./(1+N.^c(2)/c(3)))*c(4);
end

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
function perform_calibration
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
function [mlv,tpk]=tbabr_condition(fr,lv)
ts.done = 0;
spl_ref = 0.0002; % SPL reference pressure (rms Pa)
td=4/sqrt(fr); % tone-burst duration (sec)
[pa,sav,cur,ts]=tdm_init([fr;lv;td],1,4,2,ts); % initialize model
pa.ihceq=4;     % Neely synapse model
pa.hbnl=1;      % nonlinear OHC transduction
            sav=tdm_step(pa,sav,cur,ts);       % step through time
ped=sav.ped;
ped=ped - mean(ped);      % remove DC
pmn=min(ped);
pmx=max(ped);
mlv=20*log10((pmx-pmn)/(2*sqrt(2)*spl_ref));
[~,ix]=max(sav.wnr);
tpk=(ix-1)*ts.dt*ts.nw;
end % return

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

% execute ecog stimulus condition
function wnr=ecog_condition(pr,k,inv)
ts.done = 0;
lv=pr.lv(k);      % tone level (dB SPL)
fr=pr.fr(k)/1000; % tone frequency (kHz)
t1=pr.td;         % tone duration (msec)
t2=pr.ts;         % tone start    (msec)
tr=pr.tr;         % ramp duration (msec)
if (inv), t1=-t1; end % invert stimulus
stcfg=[fr;lv;t1;t2;tr];
ihcxx=pr.ihcxx(k,:);
ohcxx=pr.ohcxx(k,:);
if (length(ohcxx)>1)
    ts.ohcxx=ohcxx;
end
[pa,sav,cur,ts]=tdm_init(stcfg,1,4,2,ts); % initialize model
pa.ihceq=4;     % 4=Neely synapse model
pa.hbnl=1;      % 1=nonlinear OHC transduction
if (length(ihcxx)>1)
    if (ihcxx(2)>ihcxx(1))
        nn=1+round((pa.n-1)*(ihcxx/(10*pa.xl)));
        pa.ihcxx=false(pa.n,1);
        pa.ihcxx(nn(1):nn(2))=true;
    else
        pa.ihcxx=[];
    end
end
[sav,cur,ts,cp] = tdm_step(pa,sav,cur,ts);          % step through time
                  tdm_show(pa,sav,cur,ts,cp);       % plot saved variables
wnr=sav.wnr;
end % return

% execute loudness growth condition
function sav=ldngrw_condition(fr,lv,bw,st,rc,ew,ne)
ts.done = 0;
fr=-fr;                 % indicate narrowband stimulus
stcfg=[fr;lv;bw;st;rc]; % ldngrq stimulus configuration
[pa,sav,cur,ts]=tdm_init(stcfg,1,0,0,ts); % initialize model
pa.ihceq=4;             % Neely synapse model
pa.hbnl=1;              % nonlinear OHC transduction
if ((nargin<6)||(ew==0))
    ew=pa.ldew;
end
if ((nargin<7)||(ne==0))
    ne=pa.ldne;
end
            sav=tdm_step(pa,sav,cur,ts);       % step through time
sc=pa.ldsc;
ew=ew/(pa.xl*10);
nr=sav.nr;
nt=size(nr,1);
i1=1+round((80/120)*(nt-1)); % ramp=70, steady=50
nr=nr(i1:nt,:)';
pld = loudness(nr,ew,ne,sc);
ldso=mean(pld); % convert to sone
sav.ldso=ldso;
sav.ldcu=sone2cu(ldso);
sav.ldph=sone2phon(ldso);
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
trt=tr+50;
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
x=xx*(rms(x0)/rms(xx)).*w;
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
elseif(nf==1) % ECoG stimulus
    dur=abs(stcfg(3));
    sgn=sign(stcfg(3));
    trt=dur+20;
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
    t2=stcfg(3,k);
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
    pr=interp1(f,abs(pc_gn),fr);
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

% time-step plot
function ts=plot_ts(pa,sav,cur,ts,tt,smx)
persistent p1 p2 a2
if (ts.done), return; end % check whether ts should be plotted
if (ts.tsdsp<1), return; end
if (nargin>5)      % initialize
    ts.smx = smx / 10;
    ts.tic = tic;
    clf
    figure(ts.tsdsp)
    subplot(2,1,1)
    p1=plot(0,0,'-b');
    xlim([0 tt]);
    title(ts.lab)
    xlabel('time (ms)')
    ts.stim=pa.stim;
    switch(ts.dsp1)
        case 1, ylab='V_r';
        case 2, ylab='P_e';
        case 3, ylab='P_s';
        case 4, ylab='WNR';
    end
    xlabel('time (ms)')
    ylabel(ylab)
    subplot(2,1,2)
    p2=plot(0,0,'-b');
    a2=gca;
    if (pa.m)
        xlabel('place (mm)')
    else
        set(a2,'XtickLabel',[])
        xlabel('stapes')
    end
elseif (nargin==5) % time-step plot
    dspint=32;
    i=tt;
    if (mod(i-1,dspint)), return; end
    if (i>dspint)
        %subplot(2,1,1)
        kk=1:i;
        xx=kk*pa.ntsw*(ts.dt);
        switch(ts.dsp1)
        case 1, yy=sav.vep(kk);
        case 2, yy=sav.ped(kk);
        case 3, yy=sav.pst(kk);
        case 4, yy=sav.wnr(kk);
        end
        set(p1,'XData',xx,'YData',yy)
        %subplot(2,1,2)
        n = pa.n;
        ii = 2:(n-2); % avoid plotting endpoints
        dx = 10 * pa.xl / (n - 1);
        x = (ii-1)' * dx;
        jj = 1:pa.ldne;
        xx = (jj-1)' * (10 * pa.xl / (pa.ldne - 1));
        switch(ts.dsp2)
        case 1, y = cur.hb(ii);      % HB
        case 2, y = cur.nr(ii);      % NR
        case 3, y = cur.ld; x = xx;  % LD
        end
        set(p2,'XData',x,'YData',y)
        if (pa.m<1), ymxmn=2e-11; else, ymxmn=1e-9; end
        ymx=max(max(abs(y)),ymxmn);
        switch(ts.dsp2)
        case 1, ymn = -ymx;      % HB
        case 2, ymn = -ymx / 20; % NR
        case 3, ymn = -ymx / 20; % LD
        end
        set(a2,'YLim',[ymn ymx])
        drawnow
    end
else               % summary plot
    figure(ts.tsdsp)
    subplot(2,1,1)
    switch(ts.dsp1)
    case 1, yy=sav.vep;
    case 2, yy=sav.ped;
    case 3, yy=sav.pst;
    case 4, yy=sav.wnr;
    end
    nt=length(yy);
    xx=(1:nt)*pa.ntsw*(ts.dt);
    plot(xx,yy,'b')
    smx = max(ts.smx,max(abs(yy)));
    if (ts.dsp1>3), smn=-smx/20; else, smn=-smx; end
    if (smn==smx), smn=smx-eps; end
    axis([0 max(xx) [smn smx]*1.01]);
    subplot(2,1,2)
    kk=pa.n:-1:1;
    jj=pa.ldne:-1:1;
    switch(ts.dsp2)
    case 1, yy=sav.hb(:,kk)'; % HB
    case 2, yy=sav.nr(:,kk)'; % NR
    case 3, yy=sav.ld(:,jj)'; % LD
    end
    xc=[0 ts.ntsf*ts.dt*10];
    yc=[pa.xl*10 0];
    imagesc(xc,yc,yy)
    xlabel('time (ms)')
    if (pa.m)
        ylabel('place (mm)')
    else
        set(gca,'YtickLabel',[])
        ylabel('stapes')
    end
    if (pa.m)
        if (pa.ihcv==0), hbv='displacement'; else, hbv='velocity'; end
        switch(ts.dsp2)
        case 1, title(sprintf('hair-bundle %s',hbv)); % HB
        case 2, title('neural rate');                 % NR
        case 3, title('loudness spectrogram');        % LD
        end
    end
    if (ts.dsp2==3)
        pld = loudness(sav.nr,pa.ldew,pa.ldne,pa.ldsc);
        ldso=mean(pld);       % convert pld to sone
        ldph=sone2phon(ldso); % convert sone to phon
        fprintf('loudness = %.1f (phon)\n',ldph);
    end
end
end % return

% saved-variables plot
function plot_sv(pa,sav,ts,cp)
if (ts.done || ~ts.svdsp), return; end % check whether sv sould be plotted
fig=ts.svdsp;
if (iscell(ts.stcfg)), ts.stcfg=1; end
if (ts.stcfg), plot_tone(sav,fig);  ...          % plot tone response
else,          plot_click(pa,sav,ts,cp,fig); end % plot click response
end % return

% plot click response
function plot_click(pa,sav,ts,cp,fig)
plot_cal(pa,sav,fig)
plot_tc(pa,sav,ts,cp,fig+1)
end % return

% plot displacement tuning curves
function plot_cal(pa,sav,fig)
% Ze - middle-ear impedance (Hajicek 2024, ARO)
% Gf - middle-ear forward gain (Puria 2003, JASA)
% Zc - cochlea input impedance (Puria 2003, JASA)
ef=[0.1 0.2 0.5 1 2 5 10];                      % (kHz)
emZe=10^7*[22.8 13.3 6.5 4.8 4.1 4.0 4.9];      % (mks ohm)
emZc=1e9*[8.36 6.28 9.936 10.7 15.7 22.6 19.8]; % (mks ohm)
emGf=[-5.63 -2.77 11.0 16.6 11.7 0.25 -10.7];   % (dB)
%
f=sav.f;
f(f<1e-9)=1e-9;
nf=length(f);
ii=1:nf;
[~,~,pe,ps,vr,vs,ve]=fetch_sav(1,sav);
pc_gn=pe./vr; % ear-canal calibration
ps_gn=ps./pe; % middle-ear gain
zc=ps./(vs/pa.ast); % specific acoustic impedance
zm=pe./ve; % acoustic impedance
f=f(ii);
pc_gn=pc_gn(ii);
ps_gn=ps_gn(ii);
zc=zc(ii);
zm=zm(ii);
if (sav.nch)
    save pedcal.mat f pc_gn ps_gn
end
% plot Pe/Ve and Ps/Ps
pc_tf=dbmag(pc_gn);
ps_tf=dbmag(ps_gn);
zco=abs(zc)*10^5; % convert cgs to mks
zme=abs(zm)*10^5; % convert cgs to mks
zme(zme<1e-9)=1e-9;
figure(fig);clf
subplot(2,2,1)
loglog(f,zme,ef,emZe,'r')
xlim([0.1 20])
ylabel('Zm (mks ohm)')
title('middle-ear impedance')
subplot(2,2,2)
loglog(f,zco,ef,emZc,'r')
xlim([0.1 20])
ylabel('Zc (mks ohm)')
title('cochlea impedance')
subplot(2,2,3)
semilogx(f,ps_tf,ef,emGf,'r')
axis([0.1 20 -40 40])
xlabel('frequency (kHz)')
ylabel('Ps / Pe (dB)')
title('middle-ear gain')
subplot(2,2,4)
if (pa.m)
    f=sav.f;
    nsv=length(pa.isv);
    bf=zeros(1,nsv);
    for i=1:nsv
        [d1,d2]=fetch_sav(i,sav);
        vh=abs((d1-d2).*f);
        [~,ii]=max(log(vh));
        bf(i)=f(ii);
    end
    xl=pa.xl*10;
    xp=xl*(pa.isv-1)/(pa.n-1);
    fr=logspace(-1,log10(20),100);
    cp=f2p(fr,0.2,0.021,0.99,xl);
    semilogx(bf,xp,'b',fr,cp,'r')
    axis([0.1 20 0 xl])
    xlabel('frequency (kHz)')
    ylabel('place (mm)')
    title('frequency-place map')
else
    semilogx(f,pc_tf)
    axis([0.1 20 -30 30])
    xlabel('frequency (kHz)')
    ylabel('Pe / Vr (dB)')
    title('source calibration')
end
drawnow
end % return

% plot displacement tuning curves
function plot_tc(pa,sav,ts,cp,fig)
if (~pa.m), return; end
if (length(ts.stcfg)>1), return; end
pgd=0; % plot group delay &
f=sav.f;
s=2i*pi*f;
nf=length(f);
nsv=length(pa.isv);
spl=zeros(nf,nsv);
phv=zeros(nf,nsv);
gdv=zeros(nf,nsv);
bf=zeros(1,nsv);
qe=zeros(1,nsv);
mbf=zeros(1,nsv);
pbf=zeros(1,nsv);
dsp_ref = 140;    % dB re 1 nm displacement
spl_ref = 0.0002; % SPL reference pressure (rms Pa)
j1=ceil(nf*(0.1/max(f))); % first freqeuncy above 0.1 kHz
for i=1:nsv
    [d1,d2,pe]=fetch_sav(i,sav);
    dh = cp.hb(i,1) .* d1 + cp.hb(i,2) .* d2; % HB displacement
    ii=isnan(dh)|(abs(dh)<1e-16);
    dh(ii)=1e-16;
    vh = dh.*s;  % HB velocity (cm/s)
    pe = pe / spl_ref;
    spl(:,i)=20*log10(abs(pe./vh))+pa.hbt-dsp_ref;
    phv(:,i)=unwrap(angle(vh./pe))/(2*pi);
    gdv(:,i)=delay(vh./pe,f);
    [bfi,qei,j]=find_bf(f,-spl(:,i));
    bf(i)=bfi;
    qe(i)=qei;
    mbf(i)=spl(j,i);
    pbf(i)=phv(j,i);
    for j=j1:nf
        if (gdv(j,i)<1e-9)
            phv(j:end,i)=nan;
            gdv(j:end,i)=nan;
            break;
        end
    end
end
fprintf('BF =  ');fprintf(' %5.2f', bf );fprintf('  \n');
fprintf('Qe =  ');fprintf(' %5.2f', qe );fprintf('  \n');
fprintf('pbf=-[');fprintf(' %5.2f',-pbf);fprintf('];\n');
% recompute saved places
if ((bf(1)>0.05)&&(min(bf)==bf(1))&&(max(bf)==bf(end))&&(max(bf)~=min(bf)))
    iap=1+round(f2p(0.05,0.2,0.021,0.99,pa.n));
    isv=[iap pa.isv 1];
    fsv=[0.05 bf 20];
    isv=round(interp1(log(fsv),isv,log(0.5*2.^(-1:5))));
    fprintf('isv=[');fprintf(' %d',isv);fprintf('];\n');
end
% locate targets
fpk=500*2.^(-1:5);
pk=pk_tgt(fpk);
mpk=-pk(2,:);
% plot HB threshold tuning curves
figure(fig);clf
subplot(2,1,1)
reset_color_index;
semilogx(f,spl,fpk/1000,mpk,'ro',bf,mbf,'k.')
axis([0.1 20 0 100])
ylabel('Pe (dB SPL)')
title('HB velocity tuning curves')
text(10,15,ts.lab)
subplot(2,1,2)
reset_color_index;
if (pgd)
    %loglog(f,gdv)
    %axis([0.1 20 1 30])
    %ylabel('delay (ms)')
    %text(10,20,pa.parlab)
else
    semilogx(f,phv,bf,pbf,'k.')
    axis([0.1 20 -9 1])
    ylabel('phase (cyc)')
    text(10,-7,pa.parlab)
end
xlabel('frequency (kHz)')
drawnow
end % return

function pk=pk_tgt(flst)
fpk=500*2.^(-1:5); % 250 500 1000 2000 4000 8000 16000
xcp= [ 0.81  0.71  0.59  0.47  0.33  0.19  0.05];
mpk=-[18.20  9.70  8.80 15.00 12.30 19.10 59.10];
ppk=-[ 3.20  3.65  3.83  3.74  3.37  2.60  1.10];
xcp=interp1(fpk,xcp,flst);
mpk=interp1(fpk,mpk,flst);
ppk=interp1(fpk,ppk,flst);
% specify Zc targets
emZc=1e4*[8.36 6.28 9.936 10.7 15.7 22.6 19.8]; % (cgs ohm)
zpk=interp1(500*2.^(-1:5),emZc,flst);
pk=[xcp;mpk;ppk;zpk];
end

% plot click response
function plot_tone(sav,fig)
spl_ref = 0.0002; % SPL reference pressure (rms Pa)
f=sav.f;
nf=length(f);
ii=1:nf;
[~,~,pe]=fetch_sav(1,sav);
% plot HB threshold tuning curves
spl=dbmag(pe/spl_ref/nf/sqrt(2));
figure(fig);clf
semilogx(f,spl(ii))
axis([0.1 20 -10 90])
xlabel('frequency (kHz)')
ylabel('P_e (dB SPL)')
title('sound level at eardrum')
drawnow
end % return

function ts=plot_spec(~,~,ts)
nw = 4096;
w = hann(nw);
fs = 1 / ts.dt;
y = ts.stwav;
subplot(2,1,1)
t=linspace(0,ts.tt/1000,ts.nt);
plot(t,y);
xlim([0 max(t)])
if (ischar(ts.stcfg))
    title(ts.stcfg)
elseif (ts.stcfg==0)
    title('click')
end
subplot(2,1,2)
nov=round(nw*15/16);
S=spectrogram(y,w,nov,nw*2);
nf=512;
df=fs/nw;
S=S(1:round(nf/10),:);
t=linspace(0,ts.tt/1000,size(S,2));
f=df*(0:(nf-1));
x=f2p(f,0.2,0.021,0.99,35);
mxt=max(t);
xlb='time (ms)';
if (mxt<1)
    t=t*1000;
    mxt=max(t);
    xlb='time (ms)';
end
imagesc(t,x,abs(S))
axis([0 mxt 0 35])
xlabel(xlb)
ylabel('place (mm)')
title('Fourier spectrogram')
ts.done=1;
end

% group delay
function gd=delay(R,f)
[n,m] = size(R);
ph = unwrap(angle(R))/(2*pi);
gd = zeros([n,m]);
for k=1:m
   gd(:,k) = -cdif(ph(:,k))./cdif(f(:));
end
end % return

% centered difference
function dx=cdif(x)
n=length(x);
dx=zeros(size(x));
dx(1)=x(2)-x(1);
dx(2:(n-1))=(x(3:n)-x(1:(n-2)))/2;
dx(n)=x(n)-x(n-1);
end % return

function report_level(pa,sav,ts)
if (ts.done), return; end  % done
spl_ref = 0.0002;         % SPL reference pressure (rms Pa)
ped=sav.ped;
ped=ped - mean(ped);      % remove DC
pmn=min(ped);
pmx=max(ped);
pav=rms(ped);
lv1=20*log10(pav/spl_ref);
lv2=20*log10((pmx-pmn)/(2*sqrt(2)*spl_ref));
fprintf('eardrum pressue: %.1f (dB SPL) %.1f (dB ppeSPL) \n',lv1,lv2);
if ((pa.m)&&(pa.ihceq))
    %wnr=sav.wnr;
    %[~,ix]=max(wnr);
    %tpk=(ix-1)*ts.dt*ts.nw;
    %fprintf('WNR latency: %.1f (ms)\n',tpk);
end
end % return

%----------------------------------------------------------

function p=f2p(f,a,b,c,xl) % kHz -> mm
p = (1 - (log10((f / a) + c) / b) / 100) * xl;
end % return

%----------------------------------------------------------

function reset_color_index
if (~isoctave)
    ax=gca();
    set(ax,'ColorOrderIndex',1);
end
end % return

function o = isoctave
o = exist('OCTAVE_VERSION', 'builtin');
end % return


%==========================================================

function [d1,d2,pe,ps,vr,vs,ve,nr]=fetch_sav(i,sav)
d1 = fft(sav.d1(:,i));
d2 = fft(sav.d2(:,i));
nr = fft(sav.wnr);
vs = fft(sav.vst);
ps = fft(sav.pst);
pe = fft(sav.ped);
ve = fft(sav.ved);
vr = fft(sav.vep);
nf=length(sav.f);
ii=1:nf;
d1 = d1(ii);
d2 = d2(ii);
nr = nr(ii);
vs = vs(ii);
ps = ps(ii);
pe = pe(ii);
ve = ve(ii);
vr = vr(ii);
end % return;

function [bf,qe,mi]=find_bf(f,thr)
[mm,mi]=max(thr);
bf=f(mi);
qe=bf/trapz(f,10.^((thr-mm)/10));
end % return;

function y=dbmag(x)
y=20*log10(abs(max(x,eps)));
end % return

function x=short_chirp(nn,n1,n2,nt)
nf=nn/2+1;
mg=1;
ph=0;
X=zeros(nf,1);
for k=1:nf
    X(k)=mg*exp(-1i*ph);
    dp=(2*pi)*(n1+(k/nf)*(n2-n1))/nn;
    ph=ph+dp;
end
X=[X;zeros((nt-nn)/2,1)];
x=ffs(X);
x=x-x(1);
x=x*1.3183/rms(x); % normalize rms average
end % return

function st=flat_noise(nt,bw,cf,sr)
nf=nt/2;
ii=(0:nf)';
hb=2^(bw/2);
f0=nt*cf/sr;
f1=1+max(round(f0/hb),0);
f2=1+min(round(f0*hb),nf);
ff=(f1:f2)';
sf=zeros(size(ii));
sf(ff)=exp(2i*pi*rand(size(ff)));
st=ffs(sf);
st=st/sqrt(mean(st.^2)); % normalize to rms=1
end % return

% ffs - fast Fourier synthesize real signal
% usage: h=ffs(H)
% H - transfer function
% h - impulse response
function h=ffs(H)
m=length(H);
n=2*(m-1);
H(1,:)=real(H(1,:));
H(m,:)=real(H(m,:));
H((m+1):n,:)=conj(H((m-1):-1:2,:));
h=real(ifft(H));
end % return

function a=rms(x)
% function a=rms(x)
% x = input signal
% a = rms (root-mean-squared) amplitude of x
a=sqrt(mean(x.^2));
end % return

%==========================================================
%
function pa=modpar24(nch)
% chamber sizes
if (nch<3)
    pa=par_CEL16;
else
    pa=modpar23c3;
end
pa.m=nch;             % number of fluid chambers
pa.xp=1;              % length of excitation pattern
end

%--------------------------------------------------------------

function pa=par_CEL16
pa.parlab='CEL16';
pa.chsz=[1 1];                 % chamber size
pa.gam = 1;                    % NDR multiplier
pa.m = 1;                      % number of points across fluid
pa.n = 1401;                   % number of points along BM
pa.xl = 3.5;                   % scala length
pa.yw = 0.1;                   % scala width
pa.zh = 0.1;                   % scala height
pa.rho = 1;                    % fluid density
pa.stgain=0.0127;              % stimulus gain (target 40 dB SPL tone)
pa.isv=[1136 1005 840 655 466 273 80]; % BM locations to save
pa.hbt=8.5; pa.xtap=0.2339; pa.xtex=6; pa.mmeq=1;
pa.khe=0.0001; pa.rhe=0.0001; pa.mhe=0.01; % err=23.21 23.21
% middle-ear parameters -------------
pa.nmev=4; pa.kme=0.1; pa.rme=400; pa.mme=0.1;
pa.acp=0.5; pa.adi=0.2; pa.aed=0.5; pa.ama=0.5; pa.ast=0.03;
pa.cep=2700;  pa.rep=29.16; pa.mevgn=1.262; pa.stim=3;
pa.mcp=0.0002; pa.rcp=0.1; pa.kcp=2200; pa.rfz=0.01;
pa.mdi=0.005;  pa.rdi=100;  pa.kdi=7.7e+06;
pa.mrw=5e-3; pa.rrw=20; pa.krw=1.5e5; pa.arw=0.0625;
pa.mco=3; pa.rco=1e5;  pa.gme=0.5;
pa.kma=3e+06;  pa.rma=8000;  pa.mma=0.85;
pa.ked=7.04632e+06;  pa.red=200;  pa.med=0.02;
pa.kst=4.66022e+06;  pa.rst=8140.1;  pa.mst=0.89833;
pa.kim=3e+08;  pa.rim=28508.7; % nch=0 nme=4
% hair-cell & neural parameters --------
pa.ihceq=4; % 4 = Neely synapse model
pa.hbnl=0; pa.hbsc=0.04; pa.hbmx=6e-9; pa.ihcv=1;
pa.ihctc=0.2e-3; pa.ihcsf=1e5; pa.nrgn=2.5e6; pa.ew=0.1;
pa.ihcex=14; pa.ihcrr=17.836; pa.ihcdr=1.2637;
pa.ldew=2; pa.ldne=70; pa.ldsc=[6.913 3.823 0.454 4.499]; % 8.2%
%---------- partition impedance
pa.k1o=  2.394e+08; pa.k1e=-2.6655; pa.k1q= 0.000000; % CEL16
pa.r1o=      871.4; pa.r1e=-1.3937; pa.r1q= 0.000000;
pa.m1o=   0.007417; pa.m1e=-0.1884; pa.m1q= 0.000000;
pa.k2o=  3.036e+08; pa.k2e=-3.4762; pa.k2q= 0.000000;
pa.r2o=       1979; pa.r2e=-1.2466; pa.r2q= 0.000000;
pa.m2o=    0.03417; pa.m2e=-0.0828; pa.m2q= 0.000000;
pa.k3o=  3.151e+08; pa.k3e=-2.9092; pa.k3q= 0.000000;
pa.r3o=      1.049; pa.r3e= 0.1033; pa.r3q= 0.000000;
pa.k4o=  4.045e+08; pa.k4e=-2.7948; pa.k4q= 0.000000;
pa.r4o=          0; pa.r4e= 0.0000; pa.r4q= 0.000000;
pa.gpo=          1; pa.gpe= 0.0000; pa.gpq= 0.000000;
pa.aco=       0.01; pa.ace=-0.4000; pa.acq= 0.000000;
pa.bwo=       0.05; pa.bwe= 0.0000; pa.bwq= 0.000000;
end % return

%--------------------------------------------------------------

function pa=modpar23c3
pa.parlab='cel24c3';
pa.chsz=[0.9 0.1 1.0]; % err=93.36 93.36
pa.gam = 1;                    % NDR multiplier
pa.m = 1;                      % number of points across fluid
pa.n = 1401;                   % number of points along BM
pa.xl = 3.5;                   % scala length
pa.yw = 0.1;                   % scala width
pa.zh = 0.1;                   % scala height
pa.rho = 1;                    % fluid density
pa.stgain=0.0127;              % stimulus gain (target 40 dB SPL tone)
pa.isv=[1136 1005 840 655 466 273 80]; % BM locations to save
pa.hbt=33; pa.xtap=0.2339; pa.xtex=6; pa.mmeq=1;
pa.khe=0.0001; pa.rhe=0.0001; pa.mhe=0.01; % err=23.21 23.21
% middle-ear parameters -------------
pa.nmev=4; pa.kme=0.1; pa.rme=400; pa.mme=0.1;
pa.acp=0.5; pa.adi=0.2; pa.aed=0.5; pa.ama=0.5; pa.ast=0.03;
pa.cep=2700;  pa.rep=29.16; pa.mevgn=1.262; pa.stim=3;
pa.mcp=0.0002; pa.rcp=0.1; pa.kcp=2200; pa.rfz=0.01;
pa.mdi=0.005;  pa.rdi=100;  pa.kdi=7.7e+06;
pa.mrw=5e-3; pa.rrw=20; pa.krw=1.5e5; pa.arw=0.0625;
pa.mco=3; pa.rco=1e5;  pa.gme=0.5;
pa.kma=3e+06;  pa.rma=8000;  pa.mma=0.85;
pa.ked=7.04632e+06;  pa.red=200;  pa.med=0.02;
pa.kst=4.66022e+06;  pa.rst=8140.1;  pa.mst=0.89833;
pa.kim=3e+08;  pa.rim=28508.7; % nch=0 nme=4
% hair-cell & neural parameters --------
pa.ihceq=4; % 4 = Neely synapse model
pa.hbnl=0; pa.hbsc=0.04; pa.hbmx=6e-9; pa.ihcv=1;
pa.ihctc=0.2e-3; pa.ihcsf=1e5; pa.nrgn=2.5e6; pa.ew=0.1;
pa.ihcex=14; pa.ihcrr=17.836; pa.ihcdr=1.2637;
pa.ldew=2; pa.ldne=70; pa.ldsc=[6.913 3.823 0.454 4.499]; % 8.2%
%---------- partition impedance
pa.k1o=2.41443e+08; pa.k1e=-2.7377; pa.k1q= 0.000002;
pa.r1o=    1056.85; pa.r1e=-1.3138; pa.r1q= 0.000002;
pa.m1o=  0.0074497; pa.m1e=-0.2012; pa.m1q= 0.000002;
pa.k2o=3.05082e+08; pa.k2e=-3.5063; pa.k2q= 0.000002;
pa.r2o=    1984.21; pa.r2e=-1.2479; pa.r2q= 0.000002;
pa.m2o=  0.0360281; pa.m2e=-0.0812; pa.m2q= 0.000002;
pa.k3o=3.15079e+08; pa.k3e=-2.9231; pa.k3q= 0.000003;
pa.r3o=    1.10439; pa.r3e= 0.1065; pa.r3q= 0.000002;
pa.k4o=4.04629e+08; pa.k4e=-2.8017; pa.k4q= 0.000001;
pa.r4o=          0; pa.r4e= 0.0000; pa.r4q= 0.000000;
pa.gpo=   0.995591; pa.gpe= 0.0000; pa.gpq= 0.000003;
pa.aco= 0.00993884; pa.ace=-0.4324; pa.acq= 0.000002;
pa.bwo=  0.0517633; pa.bwe= 0.0001; pa.bwq= 0.000002;
end % return
