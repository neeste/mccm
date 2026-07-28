% ftc22.m - frequency-domain two-chamber cochlear model tuning curves
function ftc22
% compute 2DOF cochlear model 
pa=modpar22;
% chamber sizes
pa.chsz=[1 1]; % relative chamber sizes
flst=500*2.^(-1:5);
Dt=10.^(-tc_x(flst,pa.n,pa.xp,1)/20);
Dt=adj_thr(flst,Dt,1);
% compute frequency-domain cochlear model
[xx,~,~,~,Dh]=fdmod22(pa,flst);
% locate Dh peaks
bp=best_place(flst,xx*pa.xl,Dh);
% plot model results
hbt=pa.hbt;
excpat(1,xx,Dh,[-100    0],[ -7    1  ],'excitation pattern',hbt,Dt,flst);
fpmplt(2,bp,flst,pa.xl);
% compute model frequency-tuning curves
[f,epl]=tc_f(flst); % EPL average of 6 cats
epl=adj_thr(flst,epl,0);
[spl,Dh,Zc,Gf]=fdmftc(pa,f,bp);
hbdplt(3,f,Dh);
ftcplt(4,f,epl,spl);
midear(5,1000*f',pa,Zc,Gf);
end

%------------------------------------------------

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
end

%------------------------------------------------

function excpat(fig,xx,Dh,dblim,phlim,lab,hbt,Dt,flst)
tpi=2*pi;
nx=length(xx);
nf=length(flst);
db1=20*log10(abs(Dh))-hbt;
db2=20*log10(abs(Dt));
ph1=unwrap(angle(Dh))/tpi;
ph1(db1>120)=nan;
ph1=ph1-repmat(round(ph1(2,:)),nx,1);
figure(fig);clf;
subplot(211);
plot(xx,db1); hold on
ax=gca;ax.ColorOrderIndex=1;
plot(xx,db2,'--'); hold off
axis([0 1 dblim])
ylabel('magnitude (dB)')
title(lab);
subplot(212);plot(xx,ph1)
axis([0 1 phlim])
xlabel('x/L')
ylabel('phase (cyc)')
text(0.1,-6,'d1-d2')
end

%------------------------------------------------

function bp=best_place(flst,x,Dh)
% loop over frequencies
nf=length(flst);
bp=zeros(1,nf);
for k=1:nf
   [~,imx]=max(abs(Dh(:,k)));
   bp(k)=x(imx);
end
return
end

%------------------------------------------------

function fpmplt(fig,bp,flst,L)
fk=flst/1000;
cp = f2p(fk,0.1654,0.021,0.99,3.5); % man
%cp = f2p(fk,0.4560,0.021,0.80,2.5); % cat
figure(fig);clf
ax=gca;ax.ColorOrderIndex=1;
semilogx(fk,cp,'--k',fk,bp,'-o');
xlabel('frequency (kHz)')
ylabel('distance from stapes (cm)')
title('frequency-place map (man)')
axis([0.1 20 0 L])
drawnow
end

function p=f2p(fk,a,b,c,xl)
p = (1 - (log10((fk / a) + c) / b) / 100) * xl;
return
end

%------------------------------------------------

function hbdplt(fig,f,Dh)
mg=20*log10(abs(Dh))+60;
ph=angle(Dh);
ph=unwrap(ph)/(2*pi);
ph=-ph;
gd=ph;
for k=1:size(ph,2)
    gd(:,k)=-f.*cdif(ph(:,k))./cdif(f);
end
figure(fig);clf
ax=gca;ax.ColorOrderIndex=1;
ax=gca;ax.ColorOrderIndex=1;
subplot(2,1,1)
semilogx(f,mg)
axis([0.1 20 -35 55])
title('HB displacement')
ylabel('magnitude (dB)')
drawnow
subplot(2,1,2)
semilogx(f,gd)
axis([0.1 20 0 10])
xlabel('frequency (kHz)')
ylabel('phase (cyc)')
drawnow
end

function ftcplt(fig,f,ftc,spl)
figure(fig);clf
ax=gca;ax.ColorOrderIndex=1;
semilogx(f,ftc,'--'); hold on
ax=gca;ax.ColorOrderIndex=1;
semilogx(f,spl); hold off
axis([0.1 20 0 90])
title('neural threshold frequency-tuning curves (EPL,man)')
xlabel('frequency (kHz)')
ylabel('SPL (dB)')
text(0.13,5,'d1-d2')
text(10,5,'ftc22')
drawnow
end

function d=cdif(x)
d=zeros(size(x));
n=length(x);
d(1)=x(2)-x(1);
d(n)=x(n)-x(n-1);
d(2:(n-1))=x(3:n)-x(1:(n-2));
return
end
%------------------------------------------------

function [spl,Dh,Zc,Gf]=fdmftc(pa,flst,xlst)
flst=flst*1000;
nf=length(flst);
nx=length(xlst);
spl=zeros(nf,nx);
jlst=zeros(1,nx);
[xx,~,~,~,Dh,Zc,Gf]=fdmod22(pa,flst);
for k=1:nx
    [~,j]=min(abs(xlst(k)-xx(:)*pa.xl));
    jlst(k)=j;
    spl(:,k)=pa.hbt-20*log10(abs(Dh(j,:)))';
end
Dh=Dh(jlst,:)';
end

%================================================

function pa=modpar22
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
pa.xtap=0.0; pa.xtex=6;
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
end

function pa=modpar16(pa)
% ---- cel16.par ----
pa.ace=-0.4;
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
%-------- middle-ear parameters
pa.kst=5e5;
pa.rst=40;
pa.mst=0.005;
pa.ast=0.01;
pa.kma=4e5;
pa.rma=60;
pa.mma=0.01;
pa.ama=0.35;
pa.kim=4e7;
pa.rim=4e2;
pa.gmel=0.5;
pa.ked=3.1e4;
pa.red=3.2;
pa.med=0.00014;
pa.aed=0.05;
pa.kcp=2200;
pa.rcp=0.1;
pa.mcp=0.00022;
pa.acp=0.071;
pa.kdi=7.7e6;
pa.rdi=100;
pa.mdi=0.057;
pa.adi=2.5;
pa.bl=5.4e6;
pa.rvc=200;
pa.rfz=1;
%---------- partition impedance
pa.k1o=2.358e+08;
pa.k2o=3.003e+08;
pa.k3o=3.154e+08;
pa.k4o=4.112e+08;
pa.k5o=0;
pa.r1o=8.739e+02;
pa.r2o=1.928e+03;
pa.r3o=1;
pa.r4o=0;
pa.m1o=7.257e-03;
pa.m2o=3.286e-02;
pa.k1e=-2.733;
pa.k2e=-3.350;
pa.k3e=-2.942;
pa.k4e=-2.833;
pa.k5e=0;
pa.r1e=-1.305;
pa.r2e=-1.259;
pa.r3e= 0.1;
pa.r4e= 0.000;
pa.m1e=-0.061;
pa.m2e=-0.082;
pa.xtap=0.2339;
% ----
return
end

%------------------------------------------------

function [xx,Yb,Pd,Db,Dh,Zc,Gf]=fdmod22(pa,flst) % one chamber
% spatial dimensions
L=pa.xl;
N=pa.n;
dx=L/(N-1);
% physical unit reference
Pe=2.848e-4; % dB eardrum pressure
De=1e7;      % 1-nm displacement reference (cm)
nf=length(flst);
Yb=zeros(N,nf);
Pd=zeros(N,nf);
Db=zeros(N,nf);
Dh=zeros(N,nf);
Zc=zeros(1,nf);
% loop over frequencies
for k=1:nf
   f=flst(k);
   s=21*pi*f;
   f=flst(k);
   % solve for pressure
   [a,q,x,yb,hc]=mxfill(f,pa);
   pd=a\q;
   Pd(:,k)=pd/Pe; % normalize to eardrum pressure
   Yb(:,k)=yb;
   Db(:,k)=((pd.*yb)./(2i*pi*f))*De;
   Dh(:,k)=Db(:,k).*hc;
   Vs=((Pd(2,k)-Pd(1,k))/dx)/s;
   Zc(k)=Pd(1,k)/Vs;
end
xx=x/pa.xl;
Zc=Zc/pa.ast^2; % convert from mechanical to acoustic impedance
%Zc=Zc/10; % fudge factor
Gf=Pd(1,:)/2;
return
end

function [a,q,x,Yb,Hh,D]=mxfill(f,pa)
% set parameter values
H=pa.zh;
N=pa.n;
Ap=pa.bwo;
gm=pa.gam;
% useful constructs
s=2i*pi*f;
dx=pa.xl/(pa.n-1);
srd=s*pa.rho*dx;
% partition admittance
x=transpose(linspace(0,pa.xl,N));
[z1,z2,z3,z4]=imped(x,s,pa);
Hh=z2./(z2+z3);
Zm=pa.kme/s+pa.rme+pa.mme*s;
Zh=pa.khe/s+pa.rhe+pa.mhe*s;
Zb=z1+(z3-gm*z4).*Hh;
Zb(N-1)=Zh;
Yb=1./Zb;
% middle-ear boundary constraints
if (1)
    [A1,A2,A3]=midear(0,f,pa);
    A2=A2/srd;
    A1=A1-A2;
    A3=A3*1e-3;
    %A3=A3*2e-8*f*f; % cel16
else
    A1=1;
    A2=-1;
    A3=(2*srd*pa.ast./Zm)*pa.yw/2;
end
% equations for pressure
D=1;
o=ones(N,1);
%yr=(dx/H).*(Ap*Yb); ???
yr=(dx/H).*Yb;
a=[-o 2+srd*yr -o];
a(  1,1)=A1; a(1,2)=A2; % basal BC
a(N-1,N)=-1; a(N,N)=1;  % apical BC
a=spdiags(a,-1:1,N,N);
% initialize q
q=zeros(N,1);
q(1)=A3;
return
end

function [z1,z2,z3,z4]=imped(x,s,pa)
x=x.*(1+(pa.xtap*x).^pa.xtex);
q=x.^2;
k1=pa.k1o*exp(pa.k1e*x+pa.k1q*q);
r1=pa.r1o*exp(pa.r1e*x+pa.r1q*q);
m1=pa.m1o*exp(pa.m1e*x+pa.m1q*q);
k2=pa.k2o*exp(pa.k2e*x+pa.k2q*q);
r2=pa.r2o*exp(pa.r2e*x+pa.r2q*q);
m2=pa.m2o*exp(pa.m2e*x+pa.m2q*q);
k3=pa.k3o*exp(pa.k3e*x+pa.k3q*q);
r3=pa.r3o*exp(pa.r3e*x+pa.r3q*q);
k4=pa.k4o*exp(pa.k4e*x+pa.k4q*q);
r4=pa.r4o*exp(pa.r4e*x+pa.r4q*q);
z1=k1/s+r1+m1*s;
z2=k2/s+r2+m2*s;
z3=k3/s+r3;
z4=k4/s+r4;
end

function [A1,A2,A3]=midear(fig,f,pa,Zc,Gf)
fk=f/1000;
s=2i*pi*f;
dx=pa.xl/(pa.n-1);
srd=s*pa.rho*dx;
%
% middle-ear impedances
Zma=s*pa.mma+pa.rma+pa.kma;
Zim=pa.rim+pa.kim./s;
Zst=s*pa.mst+pa.rst+pa.kst./s+Zim;
Zrw=s*pa.mrw+pa.rrw+pa.krw./s+Zim;
Zme=(Zst+Zrw+Zma.*Zim./(Zma+pa.gm^2*Zim));
A1=pa.ast;
A2=Zme/srd;
A3=pa.gm*pa.aed*Zim./(Zma+pa.gm^2*Zim);
if (nargin<4) return; end
%
% cochlear input impedance
%Zco=pa.rco+pa.rrw+pa.krw./s-(pa.rco^2)./(pa.rco+pa.mco*s); % ideal
% compare with Puria, Peake, and Rosowski (1997, JASA 101,p2762)
ef1=[0.05,0.06,0.08,0.10,0.15,0.2,0.3,0.4,0.5,0.6,0.8,1.0,1.5,2,3,4,5,6,8,10];
emZco=1e4*[35,31,24,22,21,21,15,15,21,22,20,17,28,36,60,75,92,100,120,120];
%
% forward pressure gain
%Gf=A3./(pa.ast*(1-Zme./Zco));
% compare with Puria and Rosowski (1997,MOH96,p154)
ef2=[0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.5,2,3,4];
emGf=[1,3,7,13,14,14,15,16,15,14,14,13,14,13,10];
if (nargout<1)
    % plot
    figure(fig);clf
    subplot(2,1,1)
    abZco=abs(Zc);
    loglog(fk,abZco,ef1,emZco);
    xlim([0.1 20])
    title('cochlear input impedance')
    subplot(2,1,2)
    dbGf=20*log10(abs(Gf));
    semilogx(fk,dbGf,ef2,emGf);
    xlim([0.1 20])
    xlabel('frequency (Hz)')
    title('forward pressure gain')
    drawnow
end
end
