% tdm22.m - two-chamber time-domain cochlear model
%
function tdm22
tdm_init;		% set parameters and initialize state
tdm_mech;		% prepare BM & fluid equations
tdm_step;		% step through time
tdm_show;		% plot saved variables
return
%
%==========================================================
%
% set parameters and initialize state
%
function tdm_init
global pa cur sav
%
% set mechanical parameters
%pa=par_MOH90;
pa=par_CEL22;
%
% check r1/m1 ratio
lr1=log(pa.r1o) + pa.r1e + pa.r1q;
lm1=log(pa.m1o) + pa.m1e + pa.m1q;
lrm=lr1-lm1;
if (lrm>16)
    pa.nom1=1;
    fprintf('No m1: lrm=%.2f\n',lrm)
end
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
if (pa.tsdsp) figure(pa.tsdsp); clf; end
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
%
return
%-------------------------------------------
%
% prepare BM & fluid equations
%
function tdm_mech
global pa bm aa bb bst
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
b1 = -b0 ./ mf;
b2 = b0 ./ mf + [b0(2:n) ./ mf(2:n);0];
b3 = -b0 ./ mf;
% radial acceleration: zero at base
b2(1) = 0;   b3(2) = 0;
% radial acceleraton: equals longitudinal accelation
b1(n-1) = -b0(n) ./ mf(n); b2(n) = b0(n) / mf(n);
bb = spdiags([b1 b2 b3],-1:1,n,n); % tri-diagonal
% stapes BC
bst = -1 / dx;
return
%-------------------------------------------
%
% step through time
%
function tdm_step
global pa cur nxt
nxt = cur;
while (cur.cyc < pa.ncyc)
   savest(cur.stp);                         % save current state
   inctim;                                  % increment next step
  	cur_a = accel(cur, stim(cur.stp));      % get current acceleration
  	update(cur_a, cur_a);                   % integrate to next step
  	for imp=1:pa.nimp                       % improve integration?
		nxt_a = accel(nxt, stim(nxt.stp));  % get next acceleration
     	update(cur_a, nxt_a);               % integrate to next step
  	end
	cur=nxt;                                % make next state current
end
return
%-------------------------------------------
%
% plot saved variables
%
function tdm_show
global pa
if (pa.svdsp)
   nsv=length(pa.isv);
   spl=zeros(pa.ntsf/2,nsv);
   bf=zeros(1,nsv);
   ii=1:(pa.ntsf/2);
   dsp_ref = 140;  % dB re 1 nm displacement
   stm_ref = 0;    % rms stimulus level (dB SPL)
   thr_ref = dsp_ref - stm_ref; % threshold reference
   for i=1:nsv
      [f,d1,d2,pe]=fetch_sav(i);
      dh = (d1 - d2);
      thr = dbmag(dh(ii)./pe(ii)) + thr_ref;
      spl(:,i)=pa.hbt-thr;
      bf(i)=find_bf(i,f,thr);
   end
   % plot HB threshold tuning curves
   figure(pa.svdsp);clf;
   ax=gca;ax.ColorOrderIndex=1;
   semilogx(f,spl); hold on
   axis([0.1 20 0 90])
   flst=bf*1000;
   [f,spl]=tc_f(flst);
   spl=adj_thr(flst,spl,0);
   ax=gca;ax.ColorOrderIndex=1;
   semilogx(f,spl,'--'); hold off
   title('neural threshold frequency-tuning curves (EPL,man)')
   xlabel('frequency (kHz)')
   ylabel('SPL (dB)')
   text(0.13,5,'d1-d2')
   text(10,5,'tdm22')
   % estimate isv for flst
   flst=[0.25 0.5 1 2 4 8 16];
   p=polyfit(log2(bf),pa.isv,4);
   ilst=round(polyval(p,log2(flst)));
   ilst(ilst<1)=1;
   ilst(ilst>pa.n)=pa.n;
   fprintf('ilst = %d %d %d %d %d %d %d\n',ilst);
end
return
%
function spl=adj_thr(flst,spl,disp)
map=[18.2 9.7 8.8 15.0 12.3 18.1 25];
nf=length(flst);
for k=1:nf
    kk=round(log2(flst(k)/125));
    if (disp)
        mx=min(abs(20*log10(spl(:,k))));
        ad=mx-map(kk);
        spl(:,k)=spl(:,k)*10^(ad/20);
    else
        mx=min(spl(:,k));
        ad=mx-map(kk);
        spl(:,k)=spl(:,k)-ad;
    end
end
return
%
%==========================================================
%
% distribute micromechanical properties over x
%
function bm=imped(pa,bm)
%
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
%
% increment time
%
function inctim
global pa cur nxt
nxt.stp = cur.stp + 1;
nxt.cyc = cur.cyc;
if (nxt.stp >= pa.nstp)
   nxt.stp = nxt.stp - pa.nstp;
   nxt.cyc = nxt.cyc + 1;
end
return
%-------------------------------------------
% specify stimulus
%
function s=stim(i)
if (i == 0)
   s = 1;	% unit implulse
else
   s = 0;
end
return
%-------------------------------------------
% determine acceleration given state
%
function a=accel(st,stm)
global pa bm bb bst cur
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
d3 = d1 - d2;
v3 = v1 - v2;
s4 = -(bm.k4 .* d3 + bm.r4 .* v3) * pa.gam;
s3 = -(bm.k3 .* d3 + bm.r3 .* v3);
s2 = -(bm.k2 .* d2 + bm.r2 .* v2) - s3;
s1 = -(bm.k1 .* d1 + bm.r1 .* v1) + s3 - s4;
% fluid pressure
pp = press(s0, s1, stm);
% acceleration
a(i1) = bb * pp;
a(i2) = s2 ./ bm.m2;
% ME
a(ime) = (pp(2)-pp(1))*bst;
% current pressure
cur.pst = pp(1);
return
%-------------------------------------------
% integrate state to next time step
%
function update(cur_a, nxt_a)
global pa cur nxt
dt2 = pa.dt / 2;
nxt.v = cur.v + (cur_a + nxt_a) * dt2;
nxt.d = cur.d + (cur.v + nxt.v) * dt2;
return

% save data from current state
function savest(istp)
global pa cur sav
if (mod(istp,pa.ntsw) == 0)
   i = 1 + mod(round(istp / pa.ntsw), pa.ntsf);
   n = pa.n;
   sav.ped(i) = stim(istp); % assume stimulus is applied to eardrum
   sav.pst(i) = cur.pst;
   sav.d1(i,:) = cur.d(pa.isv);
   sav.d2(i,:) = cur.d(pa.isv+n);
   if (pa.tsdsp)	% time-step display
      figure(pa.tsdsp)
      ii = 1:(n-2); % avoid plotting apical end
      dx = pa.xl / (n - 1);
      x = (ii-1)' * dx;
      subplot(2,1,1)
      plot(x,cur.d(ii));
      title(sprintf('istp=%d',istp));
      drawnow;
   end
end
return
%
%==========================================================
%
% solve 1-D fluid equations
%
function p=press(s0, s1,stm)
global aa
n = length(s1);
q = s1;
stm = stm / 2;   % fudge factor ???
q(1) = stm + s0; % basal  BC: p=q(1)
q(n) = 0;        % apical BC: p(n)-p(n-1)=q(n)
p = aa \ q;
return
%
%==========================================================
%
%
% MOH90 parameters
function pa=par_MOH90
pa.gam = 1;                    % NDR multiplier
pa.m = 1;                      % number of points across fluid
pa.n = 701;                    % number of points along BM
pa.xl = 3.5;                   % scala length
pa.yw = 0.1;                   % scala width
pa.zh = 0.1;                   % scala height
pa.rho = 1;                    % fluid density
pa.bwo = 0.05;                 % BM width at base
pa.bwe = 0;                    % BM width taper
pa.isv = [555 466 378 291 205 119 31]; % BM locations to save
pa.ast= 0.01;   % area of stapes
pa.hbt=0; pa.xp=1;
pa.khe=1; pa.rhe=200; pa.mhe=1;        % parfit
pa.kme=0.01; pa.rme=1000; pa.mme=0.01; % parfit
% ---- MOH90 ----
pa.k1o=4e+08;
pa.r1o=960;
pa.m1o=0.004;
pa.k2o=2e+08;
pa.r2o=2000;
pa.m2o=0.037;
pa.k3o=2.25e+08;
pa.r3o=0;
pa.k4o=4.4e+08;
pa.r4o=0;
pa.aco=0.01;
pa.k1e=-3.0000;
pa.r1e=-1.4000;
pa.m1e=0.0000;
pa.k2e=-3.2000;
pa.r2e=-1.4000;
pa.m2e=0.0000;
pa.k3e=-3.2000;
pa.r3e=0.0000;
pa.k4e=-3.1000;
pa.r4e=0.0000;
pa.ace=0.0000;
pa.k1q=0.000000;
pa.r1q=0.000000;
pa.m1q=0.000000;
pa.k2q=0.000000;
pa.r2q=0.000000;
pa.m2q=0.000000;
pa.k3q=0.000000;
pa.r3q=0.000000;
pa.k4q=0.000000;
pa.r4q=0.000000;
pa.acq=0.000000;
return
%
% CEL22 parameters
function pa=par_CEL22
pa.gam = 1;                    % NDR multiplier
pa.m = 1;                      % number of points across fluid
pa.n = 701;                    % number of points along BM
pa.xl = 3.5;                   % scala length
pa.yw = 0.1;                   % scala width
pa.zh = 0.1;                   % scala height
pa.rho = 1;                    % fluid density
pa.bwo = 0.05;                 % BM width at base
pa.bwe = 0;                    % BM width taper
pa.isv = [562 486 408 325 235 138 34]; % BM locations to save
pa.hbt=-6; pa.xp=1;
pa.khe=0.0001; pa.rhe=0.0001; pa.mhe=0.01; % err=23.21 23.21
pa.kme=0.01; pa.rme=1000; pa.mme=0.01; % err=23.21 23.21
% middle-ear parameters
pa.mco=30; pa.rco=1.2e6; pa.rrw=2e5; pa.krw=5e7;
pa.mma=0.017; pa.rma=80; pa.kma=3e5; pa.aed=0.33;
pa.rim=400; pa.kim=5e6; pa.gm=1;
pa.mst=0.017; pa.rst=80; pa.kst=3e5; pa.ast=0.01;
pa.mrw=5e-3; pa.rrw=20; pa.krw=1.5e5; pa.arw=0.0625;
% ---- parfit ----
pa.k1o=1.054e+08;
pa.r1o=113.1;
pa.m1o=0.002;
pa.k2o=9.443e+06;
pa.r2o=15.99;
pa.m2o=0.002;
pa.k3o=1.767e+07;
pa.r3o=0.831;
pa.k4o=2.482e+07;
pa.r4o=0;
pa.aco=0.01;
pa.k1e=-0.9662;
pa.r1e=0.8793;
pa.m1e=1.7012;
pa.k2e=-2.0918;
pa.r2e=1.0379;
pa.m2e=-0.2339;
pa.k3e=-3.9166;
pa.r3e=2.5927;
pa.k4e=0.2348;
pa.r4e=0.0000;
pa.ace=0.0000;
pa.k1q=-0.290055;
pa.r1q=-0.320422;
pa.m1q=-0.583975;
pa.k2q=-0.503048;
pa.r2q=-0.540285;
pa.m2q=0.646131;
pa.k3q=0.831386;
pa.r3q=-0.359765;
pa.k4q=-0.598210;
pa.r4q=0.000000;
pa.acq=0.000000;
return
%
%==========================================================
%
function [f,d1,d2,pe,ps]=fetch_sav(i)
global pa sav
n = pa.ntsf;
f = (0:((n/2) - 1))' / (1000 * pa.dt * pa.ntsw * n);
d1 = fft(sav.d1(:,i));
d2 = fft(sav.d2(:,i));
ps = fft(sav.pst);
pe = fft(sav.ped);
return;
%
function bf=find_bf(i,f,dm)
[mm,mi]=max(dm);
bf=f(mi);
qe=bf/trapz(f,10.^((dm-mm)/10));
fprintf('isv=%d: BF=%4.2f (Hz) Qerb=%4.1f\n',i,bf,qe);
return;
%
function y=dbmag(x)
y=20*log10(abs(max(x,eps)));
return
%
