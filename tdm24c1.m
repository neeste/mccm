% tdm24c1.m - one-chamber time-domain model of cochlea
%
function tdm24c1
tdm_init(0);    % set parameters and initialize state (kHz)
tdm_mech;       % prepare BM & fluid equations
tdm_step;       % step through time
tdm_show;       % plot saved variables
return
%
%==========================================================
%
% set parameters and initialize state
%
function tdm_init(ft)
global pa cur sav stmwav
%
% set mechanical parameters
pa=par_CEL23;
%
% set timing parameters
pa.nstp = 20480;				% number of time steps per cycle
pa.ncyc = 1;					% number of cycles
pa.ntsw = 10;					% number of time steps between saves
pa.ntsf = 2048;				    % number of time steps saved
pa.dt = 2e-6;					% time step (s)
%
% set numerical integration parameter
pa.nimp = 1;					% number of "improvements"
%
% set display parameter
pa.tsdsp = 1;               % time-step display (fig. num.)
pa.svdsp = 2;               % saved variables display (fig. num.)
%
% initialize state
dof = 2;                 % degress of freedom
pa.ncp = pa.n * dof;     % number of BM variables
pa.nme = 1;              % number of ME variables
nsv = pa.ncp + pa.nme;   % number of state variables
cur.d = zeros(nsv,1);    % displacement
cur.v = zeros(nsv,1);    % velocity
cur.stp = 0;
cur.cyc = 0;
cur.ped = 0;
cur.pst = 0;
%
% initialize saved variables
nt = pa.ntsf;
nsv = length(pa.isv);
sav.ped = zeros(nt,1);
sav.pst = zeros(nt,1);
sav.dst = zeros(nt,1);
sav.d1 = zeros(nt,nsv);
sav.d2 = zeros(nt,nsv);
sav.hb = zeros(nt,pa.n);
%
% compute stimulus
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
if (pa.tsdsp)
    figure(pa.tsdsp); clf
    subplot(2,1,1)
    sav.tsxy=[0 tt min(stmwav) max(stmwav)];
    axis(sav.tsxy);
    xlabel('time (ms)')
    ylabel('P_{e}')
end
return
%-------------------------------------------
%
% prepare partition & fluid equations
%
function tdm_mech
global pa aa bb bm bc
n = pa.n;
%
% setup 1-D fluid equations
bm = imped(pa,bm);
dx = pa.xl / (n - 1);           % radial slice thickness
mf = 2 * pa.rho * dx ./ bm.ac;  % fluid inertia
mb = bm.m1 ./ (bm.bw * dx);     % BM inertia (current slice)
a1 = -mb;                       % below diagonal
a2 = mf + mb + [mb(2:n);0];     % on diagonal
a3 = -mb;                       % above diagonal
a2(1) = 1;   a3(2) = 0;         % basal BC:  p(1)=q(1)
a1(n-1) = -1; a2(n) = 1;        % apical BC: p(n)-p(n-1)=q(n)
%
aa = spdiags([a1 a2 a3],-1:1,n,n); % tri-diagonal
%
%  prepare for computationn of radial acceleration
b0 = 1 ./ (bm.bw * dx);
b0 = mb ./ bm.m1;
b1 = -b0 ./ mf;
b1 = -(mb ./ mf) ./bm.m1;
b2 = b0 ./ mf + [b0(2:n) ./ mf(2:n);0];
b3 = -b0 ./ mf;
% radial acceleration: zero at base
b2(1) = 0;   b3(2) = 0;
% radial acceleraton: equals longitudinal accelation
b1(n-1) = -b0(n) ./ mf(n); b2(n) = b0(n) / mf(n);
bb = spdiags([b1 b2 b3],-1:1,n,n); % tri-diagonal
% stapes BC
bc.st = -1 / dx;
return
%-------------------------------------------
%
% step through time
%
function tdm_step
global pa cur nxt
nxt = cur;
while (cur.cyc < pa.ncyc)
    savest(cur.stp);             % save current state
    nxt=inctim(cur,nxt,pa.nstp); % increment time to next step
    stm = stim(cur.stp);         % fetch current stimulus
  	cur_a = accel(cur,stm);      % get current acceleration
  	update(cur_a, cur_a);        % integrate to next step
  	for imp=1:pa.nimp            % improve integration?
        stm = stim(nxt.stp);     % fetch next stimulus
		nxt_a = accel(nxt,stm);  % get next acceleration
     	update(cur_a,nxt_a);     % integrate state to next step
  	end
	cur=nxt;                     % make next state current
end
return
%-------------------------------------------
%
% plot saved variables
%
function tdm_show
global pa sav
if (pa.tsdsp)
   figure(pa.tsdsp);
   subplot(2,1,2)
   kk=pa.n:-1:1;
   hb=sav.hb(:,kk)';
   xc=[0 pa.ntsf*pa.dt*pa.ntsw*1000];
   yc=[pa.xl*10 0];
   imagesc(xc,yc,hb)
   xlabel('time (ms)')
   ylabel('place (mm)')
end
if (pa.svdsp)
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
   % update saved places
   isv=[pa.n pa.isv 1];
   fsv=[0.02 bf 20];
   isv=round(interp1(fsv,isv,0.5*2.^(-1:5)));
   fprintf('isv=[');fprintf(' %d',isv);fprintf('];\n');
   % plot HB threshold tuning curves
   figure(pa.svdsp);clf;
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
   axis([0.1 20 0 xl])
end
return

function f=p2f(p,a,b,c,xl) % mm -> kHz
f = a * (10.^ (b * (1 - p / xl) * 100) - c);
return

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

% distribute micromechanical properties over x
function bm=imped(pa,bm)
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
bm.ac = pa.aco * exp(pa.ace * x + pa.acq * q);
bm.gh = pa.gpo * exp(pa.gpe * x + pa.gpq * q);
bm.bw = pa.bwo * exp(pa.bwe * x);
% basal BC
bm.kme = (pa.kme / pa.mme) * pa.yw;
bm.rme = (pa.rme / pa.mme) * pa.yw;
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

% stimulus sample
function s=stim(i)
global stmwav
s=stmwav(i+1);
return

%-------------------------------------------
% determine acceleration given state
function a=accel(st,stm)
global pa cur aa bm bb bc
a = zeros(size(st.d));
% ME pressure
ime = pa.ncp + 1;
d0 = st.d(ime);
v0 = st.v(ime);
s0 = -(bm.kme * d0 + bm.rme * v0);
% BM pressure
n = pa.n;
i1 = 1:n;
d1 = st.d(i1);
v1 = st.v(i1);
i2 = (n+1):(2*n);
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
pp = press(s0,s1,stm,aa);
% acceleration
a(i1) = bb * pp;
a(i2) = s2 ./ bm.m2;
% ME
a(ime) = (pp(2)-pp(1))*bc.st;
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
function savest(istp)
global pa cur sav
nw=pa.ntsw;
kw=mod(istp,nw);
if (kw == 0)
   i = 1 + mod(round(istp / pa.ntsw), pa.ntsf);
   n = pa.n;
   sav.ped(i) = stim(istp); % assume stimulus is applied to eardrum
   sav.pst(i) = cur.pst;
   sav.d1(i,:) = cur.d(pa.isv);
   sav.d2(i,:) = cur.d(pa.isv+n);
   sav.hb(i,:) = cur.d(1:n);
   dspint=4;
   if (pa.tsdsp && (mod(i-1,dspint)==0))	% time-step display
      figure(pa.tsdsp)
      subplot(2,1,1)
      if (i<=dspint)
          xx=[0 0];
          yy=[1 1]*sav.ped(1);
      else
          kk=(i-dspint):i;
          xx=kk*pa.ntsw*(pa.dt*1000);
          yy=sav.ped(kk);
      end
      hold on
      axis(sav.tsxy);
      plot(xx,yy,'b')
      ii = 1:(n-2); % avoid plotting apical end
      dx = 10 * pa.xl / (n - 1);
      x = (ii-1)' * dx;
      subplot(2,1,2)
      y=cur.d(ii);
      ymx=max(max(abs(y)),1e-9);
      plot(x,y);
      ylim([-ymx ymx])
      xlabel('place (mm)')
      drawnow;
   end
end
return

%==========================================================
%
% solve 1-D fluid equations
function p=press(s0, s1,stm,aa)
n = length(s1);
q = s1;
q(1) = stm + s0; % basal  BC: p=q(1)
q(n) = 0;        % apical BC: p(n)-p(n-1)=q(n)
p = aa \ q;
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
function pa=par_CEL23
pa.gam = 1;                    % NDR multiplier
pa.m = 1;                      % number of points across fluid
pa.n = 1401;                   % number of points along BM
pa.xl = 3.5;                   % scala length
pa.yw = 0.1;                   % scala width
pa.zh = 0.1;                   % scala height
pa.rho = 1;                    % fluid density
pa.bwo = 0.05;                 % BM width at base
pa.bwe = 0;                    % BM width taper
pa.isv=[1235 1054 834 655 489 306 119]; % BM locations to save
pa.hbt=-6; pa.xp=1;
pa.khe=8.814e-08; pa.rhe=1.457e-07; pa.mhe=0.05221; % err=613.7 613.6
% middle-ear parameters
pa.mco=30; pa.rco=1.2e6; pa.rrw=2e5; pa.krw=5e7;
pa.mma=0.017; pa.rma=80; pa.kma=3e5; pa.aed=0.64;
pa.rim=400; pa.kim=5e6; pa.gm=1;
pa.mst=0.017; pa.rst=80; pa.kst=3e5; pa.ast=0.03;
pa.mrw=0.005; pa.rrw=20; pa.krw=1.5e5; pa.arw=0.0625;
pa.kme=0.01529; pa.rme=355.83; pa.mme=0.012069; % err=613.6 613.1
% ---- parfit ----
pa.k1o=7.31996e+07;
pa.r1o=26.5908;
pa.m1o=0.00247877;
pa.k2o=3.48554e+07;
pa.r2o=8.27889;
pa.m2o=0.00300554;
pa.k3o=7.28254e+07;
pa.r3o=0.10505;
pa.k4o=1.55883e+07;
pa.r4o=0;
pa.aco=0.00781431;
pa.gpo=0.919832;
pa.k1e=-1.0702;
pa.r1e=0.4236;
pa.m1e=0.6666;
pa.k2e=-2.7817;
pa.r2e=2.0386;
pa.m2e=-0.1629;
pa.k3e=-5.4521;
pa.r3e=2.1064;
pa.k4e=0.0119;
pa.r4e=0.0000;
pa.ace=0.0002;
pa.gpe=-0.0013;
pa.k1q=-0.336726;
pa.r1q=-0.206915;
pa.m1q=-0.188781;
pa.k2q=-0.324471;
pa.r2q=-0.837474;
pa.m2q=0.340762;
pa.k3q=1.123875;
pa.r3q=-0.021411;
pa.k4q=-0.636990;
pa.r4q=0.000000;
pa.acq=-0.000106;
pa.gpq=-0.000608;
pa.hbt=5.3029;
return
