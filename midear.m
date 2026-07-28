function [Zma,Zim,Zst,Zrw]=midear(fig,f,pa,Zco)
s=2i*pi*f;
fk=f/1000;
%
% middle-ear impedances
Zma=s*pa.mma+pa.rma+pa.kma;
Zim=pa.rim+pa.kim./s;
Zst=s*pa.mst+pa.rst+pa.kst./s+Zim;
Zrw=s*pa.mrw+pa.rrw+pa.krw./s+Zim;
if (nargin<4) return; end
%
% cochlear input impedance
%Zco=pa.rco+pa.rrw+pa.krw./s-(pa.rco^2)./(pa.rco+pa.mco*s); % ideal
% compare with Puria, Peake, and Rosowski (1997, JASA 101,p2762)
ef1=[0.05,0.06,0.08,0.10,0.15,0.2,0.3,0.4,0.5,0.6,0.8,1.0,1.5,2,3,4,5,6,8,10];
emZco=1e4*[35,31,24,22,21,21,15,15,21,22,20,17,28,36,60,75,92,100,120,120];
%
% forward pressure gain
A1=pa.ast*(1-(Zst+Zrw+Zma.*Zim./(Zma+pa.gm^2*Zim))./Zco);
A2=pa.gm*pa.aed*Zim./(Zma+pa.gm^2*Zim);
Gf=A1./A2;
Gf=Gf/3; % fudge factor ???
% compare with Puria and Rosowski (1997,MOH96,p154)
ef2=[0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.5,2,3,4];
emGf=[1,3,7,13,14,14,15,16,15,14,14,13,14,13,10];
if (nargout<1)
    % plot
    figure(fig);clf
    subplot(2,1,1)
    abZco=abs(Zco);
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
pa.ast= 0.01;   % area of stapes
pa.hbt=0; pa.xp=1;
pa.khe=0.0001; pa.rhe=0.0001; pa.mhe=0.01; % err=23.21 23.21
% middle-ear parameters
pa.mco=30; pa.rco=1.2e6; pa.rrw=2e5; pa.krw=5e7;
pa.mma=0.017; pa.rma=80; pa.kma=3e5; pa.aed=0.5;
pa.rim=1600; pa.kim=5e6; pa.gm=2;
pa.mst=0.017; pa.rst=80; pa.kst=3e5; pa.ast=0.03;
pa.mrw=5e-3; pa.rrw=20; pa.krw=1.5e5; pa.arw=0.0625;
pa.kme=0.01; pa.rme=1000; pa.mme=0.01; % err=23.21 23.21
% ---- parfit ----
pa.k1o=1.08e+08;
pa.r1o=125.1;
pa.m1o=0.002;
pa.k2o=9.694e+06;
pa.r2o=13.66;
pa.m2o=0.002;
pa.k3o=1.839e+07;
pa.r3o=1.022;
pa.k4o=2.328e+07;
pa.r4o=0;
pa.aco=0.01;
pa.k1e=-0.9921;
pa.r1e=0.9516;
pa.m1e=1.7687;
pa.k2e=-2.2658;
pa.r2e=1.1247;
pa.m2e=-0.2391;
pa.k3e=-3.6825;
pa.r3e=2.5807;
pa.k4e=0.1818;
pa.r4e=0.0000;
pa.ace=0.0000;
pa.k1q=-0.279247;
pa.r1q=-0.361218;
pa.m1q=-0.494342;
pa.k2q=-0.457552;
pa.r2q=-0.612309;
pa.m2q=0.652367;
pa.k3q=0.764653;
pa.r3q=-0.385417;
pa.k4q=-0.577444;
pa.r4q=0.000000;
pa.acq=0.000000;
return
end

%------------------------------------------------

function [xx,Yb,Pd,Db,Dh,Zc]=fdmod22(pa,flst) % one chamber
% spatial dimensions
L=pa.xl;
N=pa.n;
dx=L/(N-1);
% physical unit reference
Pe=2.848e-4; % 0 dB SPL at the eardrum
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
% equations for pressure
D=1;
o=ones(N,1);
%yr=(dx/H).*(Ap*Yb); ???
yr=(dx/H).*Yb;
a=[-o 2+srd*yr -o];
a(  1,1)=1;  a(1,2)=-1; % basal BC
a(N-1,N)=-1; a(N,N)=1;  % apical BC
a=spdiags(a,-1:1,N,N);
% initialize q
q=zeros(N,1);
q(1)=(2*srd*pa.ast./Zm)*pa.yw*sqrt(2);
return
end

function [z1,z2,z3,z4]=imped(x,s,pa)
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

