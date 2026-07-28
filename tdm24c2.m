% tdm24c2.m - two-chamber time-domain model of cochlea
%
function tdm24c2
tdm_init(1,1);  % initialize model parameters & state (kHz,nch)
tdm_mech;       % prepare BM & fluid mechanical equations
tdm_step;       % step through time
tdm_show;       % plot saved variables
return
%==========================================================
%
% configure options, initialize model parameters and state
%
function tdm_init(ft,nch)
global pa cur sav
%
% specify mechanical parameters
pa=modpar24(nch);
fprintf('nch = %d\n',nch)
%
% specify timing parameters
pa.nstp = 20480;         % number of time steps per cycle
pa.ncyc = 1;             % number of cycles
pa.ntsw = 10;            % number of time steps between saves
pa.ntsf = 2048;          % number of time steps saved
pa.dt = 2e-6;            % time step (s)
%
% set numerical integration parameter
pa.nimp = 1;             % number of "improvements"
%
% initialize state
dof = 2;                 % degress of freedom
pa.ncp = pa.n * dof;     % number of BM variables
nsv = pa.ncp + pa.nmev;  % number of state variables
cur.d = zeros(nsv,1);    % displacement
cur.v = zeros(nsv,1);    % velocity
cur.stp = 0;
cur.cyc = 0;
cur.wnr = 0;
cur.pst = 0;
cur.ped = 0;
cur.vep = 0;
cur.stm = 0;
cur.qme = zeros(pa.nmev,1);
%
% initialize saved variables
nt = pa.ntsf;
nsv = length(pa.isv);
sav.vep = zeros(nt,1);
sav.ped = zeros(nt,1);
sav.pst = zeros(nt,1);
sav.dst = zeros(nt,1);
sav.d1 = zeros(nt,nsv);
sav.d2 = zeros(nt,nsv);
sav.hb = zeros(nt,pa.n);
%
% initialize stimulus & graphics
ts_init(1,2,ft);  % figure numbers, total time
return
%
%==========================================================
%
% prepare mechanics, set up equations of motion
%
function tdm_mech
midear  % middle-ear and coupler
cochlea % cochlear partition & fluid
neural  % neural transduction
return

function midear
global pa me
if (pa.nmev==1) % simplified middle-ear
    me.smek = (pa.kme / pa.mme) * pa.yw;
    me.smer = (pa.rme / pa.mme) * pa.yw;
    me.smes = 1;
    if (pa.stim==3) 
        pa.stim=2; 
        fprintf('moved earphone stimulus to eardrum\n')
    end
    me.ame = 0;
    me.pedm = 0;
    me.pedr = 0;
    me.pedk = 0;
    me.smem = (pa.n - 1) / pa.xl;
    return
end
nme=pa.nmev; % number of middle-ear components
% middle-ear-component indices
iep=1; % earphone diaphragm
ied=2; % eardrum (non-propagating motion)
ima=3; % malleus
ist=4; %s stapes
% middle-ear area
me.ame(iep)=pa.adi;
me.ame(ied)=pa.acp;
me.ame(ima)=pa.ama;
me.ame(ist)=pa.ast;
% middle-ear useful constructs
cep = 0.1 * (pa.bl / pa.rvc);
rep = 1e-9 * ((pa.bl * pa.bl) / pa.rvc);
ace = pa.acp / pa.aed;
adc = pa.adi / pa.acp;
aec = pa.aed / pa.acp;
amc = pa.ama / pa.acp;
amx = pa.ama / pa.aed;
gme = pa.gmel;
ma = pa.mcp;
ra = pa.rcp + pa.rfz;
rb = pa.rcp * adc;
ka = pa.kcp;
kb = pa.kcp * adc;
% set all middle-ear vaiables to zero
me.pedm=zeros(nme,1); me.pedr=zeros(nme,1); me.pedk=zeros(nme,1); 
me.smes=zeros(nme,1);
me.smem=zeros(nme,nme); me.smer=zeros(nme,nme); me.smek=zeros(nme,nme);
% set middle-ear mass
me.smem(iep,iep) = pa.mdi;
me.smem(ied,ied) = ace * pa.med + aec * ma;
me.smem(ied,ima)= amx * pa.med;
me.smem(ima,ied) = amc * ma;
me.smem(ima,ima) = pa.mma;
me.smem(ist,ist) = pa.mst;
% set middle-ear damping
me.smer(iep,iep) = pa.rdi + rep + pa.rcp * adc^2;
me.smer(iep,ied) = -pa.rcp * adc;
me.smer(ied,iep) = -aec * rb;
me.smer(ied,ied) = ace * pa.red + aec * ra;
me.smer(ied,ima) = -amx * pa.red;
me.smer(ima,iep) = -amc * rb;
me.smer(ima,ied) = amc * ra;
me.smer(ima,ima) = pa.rma +  gme^2 * pa.rim;
me.smer(ima,ist) = -gme * pa.rim;
me.smer(ist,ima) = -gme * pa.rim;
me.smer(ist,ist) = pa.rst + pa.rim;
% set middle-ear stiffness
me.smek(iep,iep) = pa.kdi + adc^2 * pa.kcp;
me.smek(iep,ied) = -adc * pa.kcp;
me.smek(ied,iep) = -aec * kb;
me.smek(ied,ied) = ace * pa.ked + aec * ka;
me.smek(ied,ima) = -amx * pa.ked;
me.smek(ima,iep) = -amc * kb;
me.smek(ima,ied) = amc * ka;
me.smek(ima,ima) = pa.kma +  gme^2 * pa.kim;
me.smek(ima,ist) = -(gme * pa.kim);
me.smek(ist,ima) = -(gme * pa.kim);
me.smek(ist,ist) = pa.kst + pa.kim;
% divde sme by ame
me.smem=me.smem./me.ame;
me.smer=me.smer./me.ame;
me.smek=me.smek./me.ame;
% invert mass matrix (smem)
a = me.smem(ied,ied);
b = me.smem(ied,ima);
c = me.smem(ima,ied);
d = me.smem(ima,ima);
det = a * d - b * c;
if (det < 1e-9) % ill-conditioned ?
    error('bad middle-ear paramters');
end
me.smem(iep,iep) = 1  / me.smem(iep,iep);
me.smem(ied,ied) = d / det;
me.smem(ied,ima) = -b / det;
me.smem(ima,ied) = -c / det;
me.smem(ima,ima) = a / det;
me.smem(ist,ist) = 1  / me.smem(ist,ist);
% apply stimulus
if (pa.stim==3)     % --> earphone
    me.smes(iep) = (cep / pa.adi);
elseif (pa.stim==2) % -->  eardrum
    me.smes(ied) = 1;
    me.smes(ima) = 1;
end
me.smes = me.smes * pa.mestgn;
me.peds = 0;
me.pedm(ied) = (pa.med / pa.aed) * (pa.acp / pa.aed);
me.pedm(ima) = -(pa.med / pa.aed) * (pa.ama / pa.aed);
me.pedr(ied) = (pa.red / pa.aed) * (pa.acp / pa.aed);
me.pedr(ima) = -(pa.red / pa.aed) * (pa.ama / pa.aed);
me.pedk(ied) = (pa.ked / pa.aed) * (pa.acp / pa.aed);
me.pedk(ima) = -(pa.ked / pa.aed) * (pa.ama / pa.aed);
% stapes boundary condition
dx = pa.xl / (pa.n - 1);
astom = me.ame(ist) * me.smem(ist,ist);
me.alfx = ((4 * pa.rho * dx) / pa.aco) * astom;
return

function cochlea
global pa aa bb bm
n = pa.n;
m = pa.m;
nm = n * m;
%
% setup m-chamber fluid equations
a1=zeros(nm,1);a2=zeros(nm,1);a3=zeros(nm,1);
bm = imped(pa);
kk = 2:(n-1);
dx = pa.xl / (n - 1);                % radial slice thickness
ar = bm.bw * dx * 2;                 % BM area
ap = bm.m1 .* bm.ac / (pa.rho * dx); % BM/fluid mass-ratio
a2(1) = 1;                           % basal BC:  p(1)=q(1)
a1(kk) = -ap(kk);                    % left of diagonal
a3(kk) = -ap(kk+1);                  % right of diagonal
a2(kk) = ar(kk)+ap(kk)+ap(kk+1);     % diagonal element
a1(n) = -1;                          % apical BC: p(n)-p(n-1)=q(n)
a2(n) = 1;
%
% format matrix A
aa=xpnd_a(a1,a2,a3,m,n);
%
%  prepare for computation of radial acceleration
ra = -bm.ac./((2*pa.rho*dx)*(bm.bw*dx));
% format matrix B
bb=xpnd_b(ra,m,n);
return

function neural
return

%-------------------------------------------

function aa=xpnd_a(a1,a2,a3,m,n)
nm=n*m;
a1x=zeros(nm,1);
a2x=zeros(nm,1);
a3x=zeros(nm,1);
for k=1:n
    kk = (k - 1) * m;
    a1x(kk+1)=a1(k);
    a2x(kk+1)=a2(k);
    a3x(kk+1)=a3(k);
    for j=2:m
        a2x(kk+j)=1;
    end
end
% format as sparse 3m-diagonal matrix
a1 = [a1x((m+1):nm);zeros(m,1)];
a2 = a2x;
a3 = [zeros(m,1);a3x(1:(nm-m))];
aa = spdiags([a1 a2 a3],[-m 0 m],nm,nm);
if (m>1)
    for k=2:n
        for j=2:m
            kk = k*m;
            aa(kk,kk-(j-1))=1;
        end
    end
end
return

function bb=xpnd_b(ra,m,n)
kk=2:(n-1);
b1=zeros(size(ra));b3=zeros(size(ra));
b1(kk) = ra(kk-1);
b3(kk) = ra(kk+1);
b2 = -(b1+b3);
b2(1) = -b3(1); % base
b2(n) = -b1(1); % apex
nm=n*m;
b1x=zeros(nm,1);
b2x=zeros(nm,1);
b3x=zeros(nm,1);
for k=1:n
    kk = (k - 1) * m;
    b1x(kk+1)=b1(k);
    b2x(kk+1)=b2(k);
    b3x(kk+1)=b3(k);
end
% format as sparse 3m-diagonal matrix
ra = [b1x((m+1):nm);zeros(m,1)];
b2 = b2x;
b3 = [zeros(m,1);b3x(1:(nm-m))];
bb = spdiags([ra b2 b3],[-m 0 m],nm,nm);
return

function qx=xpnd_q(qst,s1,m,n)
global bm
qx=zeros(n*m,1);
q = s1 .* (bm.bw .* bm.ac);
q(1) = qst; % basal  BC: p=q(1)
q(n) = 0;   % apical BC: p(n)-p(n-1)=q(n)
if (m==1)
    qx = q;
elseif (m==2)
    i2 = (1:n)*2;
    i1 = i2-1;
    qx(i1) = q;
    qx(i2) = q;
end
return

function p=fold_p(px,m,n)
if (m==1)
    p=px;
elseif (m==2)
    i2=(1:n)*2;
    i1=i2-1;
    p=px(i1)-px(i2); % (ST-SV) pressure difference 
end
return

%-------------------------------------------
%
% step through time
%
function tdm_step
global pa cur nxt
nxt = applst(cur);                 % apply stimulus to initial state
while (cur.cyc < pa.ncyc)
    savest;                        % save current state
    nxt = inctim(cur,nxt,pa.nstp); % increment time to next step
  	cur_a = accel(cur);            % get current acceleration
  	update(cur_a, cur_a);          % integrate to next step
    nxt = applst(nxt);             % apply stimulus to next state
  	for imp=1:pa.nimp              % improve integration?
		nxt_a = accel(nxt);        % get next acceleration
     	update(cur_a,nxt_a);       % integrate state to next step
  	end
	cur=nxt;                       % make next state current
end
return
%-------------------------------------------
%
% plot saved variables
%
function tdm_show
plot_ts; % time-step plot (summary)
plot_sv; % saved-variables plot
return

%==========================================================

% distribute micromechanical properties over x
function bm=imped(pa)
n = pa.n;
dx = pa.xl / (n - 1);
x = (0:(n - 1))' * dx;
q = x .^ 2;
bm.k1 = pa.k1o * exp(pa.k1e * x + pa.k1q * q);
bm.r1 = pa.r1o * exp(pa.r1e * x + pa.r1q * q);
bm.m1 = pa.m1o * exp(pa.m1e * x + pa.m1q * q);
bm.k2 = pa.k2o * exp(pa.k2e * x + pa.k2q * q);
bm.r2 = pa.r2o * exp(pa.r2e * x + pa.r2q * q);
bm.m2 = pa.m2o * exp(pa.m2e * x + pa.m2q * q);
bm.k3 = pa.k3o * exp(pa.k3e * x + pa.k3q * q);
bm.r3 = pa.r3o * exp(pa.r3e * x + pa.r3q * q);
bm.k4 = pa.k4o * exp(pa.k4e * x + pa.k4q * q);
bm.r4 = pa.r4o * exp(pa.r4e * x + pa.r4q * q);
bm.gh = pa.gpo * exp(pa.gpe * x + pa.gpq * q);
bm.ac = pa.aco * exp(pa.ace * x + pa.acq * q);
bm.bw = pa.bwo * exp(pa.bwe * x + pa.bwq * q);
% apical BC
ihe=n;
bm.k1(ihe) = pa.khe;
bm.r1(ihe) = pa.rhe;
bm.m1(ihe) = pa.mhe;
bm.k2(ihe) = 0;
bm.r2(ihe) = 0;
bm.m2(ihe) = 0;
bm.k3(ihe) = 0;
bm.r3(ihe) = 0;
bm.k4(ihe) = 0;
bm.r4(ihe) = 0;
return

% increment time
function nxt=inctim(cur,nxt,nstp)
nxt.stp = cur.stp + 1;
nxt.cyc = cur.cyc;
if (nxt.stp >= nstp)
   nxt.stp = nxt.stp - nstp;
   nxt.cyc = nxt.cyc + 1;
end
return
%-------------------------------------------

% apply stimulus to given state
function st=applst(st)
global pa stmwav
stm = stmwav(st.stp+1);
switch (pa.stim)
    case 1, st.pst=stm;
    case 2, st.pst=stm;
    case 3, st.vep=stm;
end
st.stm=stm;
st.qme(pa.stim)=stm;
return

%-------------------------------------------
% determine acceleration given state
function a=accel(st)
global pa cur aa bm bb me
a = zeros(size(st.d));
stm = st.stm;
% ME pressure
nme = pa.nmev;
ime = pa.ncp + (1:nme);
d0 = st.d(ime);
v0 = st.v(ime);
s0 = (me.smek * d0 + me.smer * v0) - stm * me.smes;
ist = nme;  % stapes index
% BM pressure
m = pa.m;
n = pa.n;
i1 = 1:n;
i2 = (n+1):(2*n);
d1 = st.d(i1);
v1 = st.v(i1);
d2 = st.d(i2);
v2 = st.v(i2);
% mmeq = 1
d3 = bm.gh .* d1 - d2;
v3 = bm.gh .* v1 - v2;
s4 = -(bm.k4 .* d3 + bm.r4 .* v3) * pa.gam;
s3 = -(bm.k3 .* d3 + bm.r3 .* v3);
s2 = -(bm.k2 .* d2 + bm.r2 .* v2) - s3;
s1 = -(bm.k1 .* d1 + bm.r1 .* v1) + s3 - s4;
% fluid pressure
if (pa.stim==1) % stapes pressure stimulus
    qst = stm;
else           % compute force on stapes
    prw = (pa.rrw * v0(ist) + pa.krw * d0(ist)) * me.ame(ist);
    qst = s0(ist) + prw;
end
if (pa.n>1)
    pp = press(qst,s1,aa,m,n); % whole cochlea
else
    rco = pa.rco * me.ame(ist);
    pp = qst + rco * v0(ist);  % cochlear input impedance only
end
% ME acceleration
if (nme==1)
    %s0 = -pp(1+m);
end
a0 = -me.smem * (pp(1) + s0(ist) + prw);
% update eardrum pressure, if not specified
if (pa.stim~=2)
    cur.ped = me.pedm .* a0 + me.pedr .* v0 + me.pedk .* d0;
end
a(ime)=a0; % return middle-ear acceleration
% cochlea acceleration
bbpp = fold_p(bb * pp,m,n);
a(i1) = bbpp;
a(i2) = s2 ./ bm.m2;
% current pressure
cur.pst = pp(1);
return
%-------------------------------------------
% integrate state to next time step
function update(cur_a, nxt_a)
global pa cur nxt
dt2 = pa.dt / 2;
nxt.v = cur.v + (cur_a + nxt_a) * dt2;
nxt.d = cur.d + (cur.v + nxt.v) * dt2;
return

% save data from current state
function savest
global pa cur sav
istp = cur.stp;
nw=pa.ntsw;
kw=mod(istp,nw);
if (kw == 0)
   i = 1 + mod(round(istp / pa.ntsw), pa.ntsf);
   n = pa.n;
   sav.vep(i) = cur.vep;
   sav.ped(i) = cur.stm;
   sav.pst(i) = cur.pst;
   sav.d1(i,:) = cur.d(pa.isv);
   sav.d2(i,:) = cur.d(pa.isv+n);
   sav.hb(i,:) = cur.d(1:n);
   % time-step display
   plot_ts(i); 
end
return

%==========================================================
%
% solve fluid equations
function p=press(qst,s1,aa,m,n)
q=xpnd_q(qst,s1,m,n);
p = aa \ q;
return

%==========================================================

% initialize time-step + saved-state displays
function ts_init(tsfig,svfig,ft)
global pa
pa.tsdsp = tsfig;           % time-step figure number
pa.svdsp = svfig;           % saved-variables figure number
pa.styp = ft;               % stimulus type
[tt,smn,smx]=stim_init(ft); % initialize stimulus
plot_ts(tt,smn,smx);        % initialize time-step plot
return

function [tt,smn,smx]=stim_init(ft)
global pa stmwav
dt=2e-5; % 20 usec
nn=round(pa.nstp*pa.dt/dt);
tt=nn*dt*1000;
nt=pa.nstp;
if (ft) % tone stimulus
    yy=sin(2*pi*ft*linspace(0,tt,pa.nstp));
else      % chirp chirp stimulus
    n1=1+round(0.001/dt);
    n2=1+round(0.002/dt);
    yy=short_chirp(nn,n1,n2,nt);
end
stmwav=yy/rms(yy)/sqrt(length(yy));
smn=min(stmwav);
smx=max(stmwav);
return

% time-step plot
function plot_ts(tt,smn,smx)
global pa sav cur
if (~pa.tsdsp) return; end
figure(pa.tsdsp)
if (nargin>1)      % initialize
    clf
    subplot(2,1,1)
    sav.tsxy=[0 tt smn smx];
    axis(sav.tsxy);
    title(sprintf('nch=%d',pa.m))
    xlabel('time (ms)')
    switch(pa.stim)
        case 1, ylab='P_s';
        case 2, ylab='P_e';
        case 3, ylab='V_e';
    end
    xlabel('time (ms)')
    ylabel(ylab)
elseif (nargin==1) % time-step plot
    dspint=4;
    i=tt;
    if (mod(i-1,dspint)) return; end
    subplot(2,1,1)
    if (i>dspint)
        kk=(i-dspint):i;
        switch(pa.stim)
        case 1, dsp=sav.pst(kk);
        case 2, dsp=sav.ped(kk);
        case 3, dsp=sav.vep(kk);
        end
        xx=kk*pa.ntsw*(pa.dt*1000);
        yy=dsp;
        hold on
        plot(xx,yy,'b')
    end
    subplot(2,1,2)
    n = pa.n;
    ii = 2:(n-2); % avoid plotting endpoints
    dx = 10 * pa.xl / (n - 1);
    x = (ii-1)' * dx;
    y=cur.d(ii);
    ymx=max(max(abs(y)),1e-9);
    plot(x,y);
    ylim([-ymx ymx])
    xlabel('place (mm)')
    drawnow;
else               % summary plot
    subplot(2,1,2)
    kk=pa.n:-1:1;
    hb=sav.hb(:,kk)';
    xc=[0 pa.ntsf*pa.dt*pa.ntsw*1000];
    yc=[pa.xl*10 0];
    imagesc(xc,yc,hb)
    xlabel('time (ms)')
    ylabel('place (mm)')
end
return

% saved-variables plot
function plot_sv
global pa
if (~pa.svdsp) return; end       % check figure number
fig=pa.svdsp;
if (~pa.styp) plot_click(fig);...  % plot click response
else          plot_tone(fig); end; % plot tone response
return

% plot click response
function plot_click(fig)
plot_tc(fig)
return

% plot displacement tuning curves
function plot_tc(fig)
global pa sav
nsv=length(pa.isv);
spl=zeros(pa.ntsf/2,nsv);
bf=zeros(1,nsv);
qe=zeros(1,nsv);
ii=1:(pa.ntsf/2);
dsp_ref = 140;  % dB re 1 nm displacement
stm_ref =  30;  % rms stimulus level (dB SPL)
thr_ref = dsp_ref - stm_ref; % threshold reference
for i=1:nsv
    [f,d1,d2,pe]=fetch_sav(i,pa,sav);
    dh = (d1 - d2);
    thr = dbmag(dh(ii)./pe(ii)) + thr_ref;
    spl(:,i)=pa.hbt-thr;
    [bfi,qei]=find_bf(f,thr);
    bf(i)=bfi;
    qe(i)=qei;
end
fprintf('BF =');fprintf(' %5.2f',bf);fprintf('\n');
fprintf('Qe =');fprintf(' %5.2f',qe);fprintf('\n');
% recompute saved places
f1=0.05;
i1=1+round(f2p(f1,0.2,0.021,0.99,pa.n));
isv=[i1 pa.isv 1];
fsv=[f1 bf 20];
isv=round(interp1(log(fsv),isv,log(0.5*2.^(-1:5))));
fprintf('isv=[');fprintf(' %d',isv);fprintf('];\n');
% plot HB threshold tuning curves
figure(fig);clf
subplot(2,1,1)
reset_color_index;
semilogx(f,spl); hold on
flst=bf*1000;
if (exist('tc_f.m','file'))
    [f,spl]=tc_f(flst);
    spl=adj_thr_f(flst,spl,0);
    reset_color_index;
    semilogx(f,spl,'--')
    title('neural threshold frequency-tuning curves (EPL,man)')
else
    title('HB iso-displacement frequency-tuning curves')
end
hold off
axis([0.1 20 0 90])
ylabel('SPL (dB)')
text(10,5,'cel23')
% estimate isv for flst
subplot(2,1,2)
xl=pa.xl*10;
xmod=xl*(pa.isv-1)/(pa.n-1);
fmod=bf;
xman=linspace(0,xl,256);
fman=p2f(xman,0.2,0.021,0.99,xl); % man
semilogx(fmod,xmod,'ro',fman,xman,'b')
xlabel('frequency (kHz)')
ylabel('place (mm)')
legend('model','man')
text(0.2,xl/10,sprintf('nch=%d',pa.m))
axis([0.1 20 0 xl])
return

% plot click response
function plot_tone(fig)
return

%----------------------------------------------------------

function reset_color_index
if (~isoctave)
    ax=gca();
    set(ax,'ColorOrderIndex',1);
end
return

function o = isoctave
o = exist('OCTAVE_VERSION', 'builtin');
return

function spl=adj_thr_f(flst,spl,disp)
map=[18.2 9.7 8.8 15.0 12.3 19.1 59.1];
nf=length(flst);
for k=1:nf
    kk=round(log2(flst(k)/125));
    if (disp)
        mx=min(abs(20*log10(spl(:,k))));
        ad=mx-map(kk);
        spl(:,k)=spl(:,k)*10^(ad/20);
    elseif (isnan(kk)||(kk<1)||(kk>nf))
        spl(:,k)=0;
    else
        mx=min(spl(:,k));
        ad=mx-map(kk);
        spl(:,k)=spl(:,k)-ad;
    end
end
return

%==========================================================

function [f,d1,d2,pe,ps]=fetch_sav(i,pa,sav)
n = pa.ntsf;
f = (0:((n/2) - 1))' / (1000 * pa.dt * pa.ntsw * n);
d1 = fft(sav.d1(:,i));
d2 = fft(sav.d2(:,i));
ps = fft(sav.pst);
pe = fft(sav.ped);
return;

function [bf,qe]=find_bf(f,thr)
[mm,mi]=max(thr);
bf=f(mi);
qe=bf/trapz(f,10.^((thr-mm)/10));
return;

function y=dbmag(x)
y=20*log10(abs(max(x,eps)));
return

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
return

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
return

function a=rms(x)
% function a=rms(x)
% x = input signal
% a = rms (root-mean-squared) amplitude of x
a=sqrt(mean(x.^2));
return

%==========================================================
%
function pa=modpar24(nch)
pa=modpar24c2;
pa.m = nch;
return

function pa=modpar24c2
pa.chsz=[1 1];                 % two-chamber sizes
% CEL23 parameters --------------------------------------------
pa.gam = 1;                    % NDR multiplier
pa.n = 1401;                   % number of points along BM
pa.xl = 3.5;                   % scala length
pa.yw = 0.1;                   % scala width
pa.zh = 0.1;                   % scala height
pa.rho = 1;                    % fluid density
pa.bwo = 0.05;                 % BM width at base
pa.bwe = 0;                    % BM width taper
pa.isv=[1185 1046 838 657 489 307 121]; % BM locations to save
pa.hbt=-6; pa.xp=1;
pa.khe=0.0001; pa.rhe=0.0001; pa.mhe=0.01; % err=23.21 23.21
% middle-ear parameters
pa.nmev = 1; % middle-ear components [4=complete]
pa.acp=0.5; pa.adi=0.2; pa.aed=0.33; pa.ama=0.35; pa.ast=0.01;
pa.bl=5.4e+06; pa.rvc=200; pa.mestgn=1; pa.stim=3;
pa.mcp=0.0002; pa.rcp=0.1; pa.kcp=2200; pa.rfz=1;
pa.mdi=0.005; pa.rdi=100; pa.kdi=7.7e+06;
pa.med=0.0001; pa.red=3.2; pa.ked=3.1e+04;
pa.mma=0.017; pa.rma=80; pa.kma=3e5;
pa.rim=400; pa.kim=5e6; pa.gmel=0.5;
pa.mst=0.017; pa.rst=80; pa.kst=3e5;
pa.mco=30; pa.rco=1.2e6; pa.rrw=2e5; pa.krw=5e7;
pa.mrw=5e-3; pa.rrw=20; pa.krw=1.5e5; pa.arw=0.0625;
pa.kme=0.01476; pa.rme=400.27; pa.mme=0.010887; % err=184.3 184.3
% ---- parfit ----
pa.k1o=7.31986e+07;
pa.r1o=32.3202;
pa.m1o=0.00214215;
pa.k2o=3.45465e+07;
pa.r2o=9.94087;
pa.m2o=0.00297927;
pa.k3o=7.30666e+07;
pa.r3o=0.127682;
pa.k4o=1.70562e+07;
pa.r4o=0;
pa.gpo=0.987832;
pa.aco=0.00992988;
pa.bwo=0.0495006;
pa.k1e=-1.0888;
pa.r1e=0.4234;
pa.m1e=0.6531;
pa.k2e=-2.8978;
pa.r2e=2.2038;
pa.m2e=-0.1484;
pa.k3e=-5.4991;
pa.r3e=1.9980;
pa.k4e=0.0148;
pa.r4e=0.0000;
pa.gpe=-0.0012;
pa.ace=-0.1000;
pa.bwe=0.1000;
pa.k1q=-0.328961;
pa.r1q=-0.251789;
pa.m1q=-0.212382;
pa.k2q=-0.195127;
pa.r2q=-0.975343;
pa.m2q=0.323923;
pa.k3q=1.096276;
pa.r3q=-0.027639;
pa.k4q=-0.656912;
pa.r4q=0.000000;
pa.gpq=-0.000812;
pa.acq=0.000000;
pa.bwq=0.000000;
pa.hbt=6.0362;
return
