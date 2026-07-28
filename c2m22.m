% c2m22 2-chamber cochlear model - 2022
function c2m22
fn='c2m22';nfig=0;
% fetfch model parameters
pa=modpar22;
pa.chsz=[1 1]; % relative chamber sizes
flst=400*(sqrt(2).^(0:9));
% compute frequency-domain cochlear model
[xx,Yb,Pd,Db,Dh]=fdmod22(pa,flst);
% plot model results
fdm_plt(flst,xx,Yb,Pd,Db,Dh,pa,nfig);
return
end

function fdm_plt(flst,xx,Yb,Pd,Db,Dh,pa,nfig)
L=pa.xl;
x=xx*L;
% prepare figures
mpfrm(1,[0 L -100   0],[0 L -0.3 0.3],'BM admittance');
mpfrm(2,[0 L  -80  20],[0 L -9.0 1.0],'BM presure difference');
mpfrm(3,[0 L -110 -10],[0 L -9.0 1.0],'BM displacement re stapes');
mpfrm(4,[0 L -110 -10],[0 L -9.0 1.0],'HB displacement re stapes');
mpfrm(5,[0 L  -15  15],[0 L -0.3 0.3],'HB transfer');
% loop over frequencies
nf=length(flst);
cp=zeros(1,nf);
for k=1:nf
   f=flst(k);
   pd=Pd(:,k);
   yb=Yb(:,k);
   d1=Db(:,k);
   d2=Dh(:,k);
   hh=d2./d1;
   % plot data
   mpplt(1,x,yb);
   mpplt(2,x,pd);
   mpplt(3,x,d1);
   mpplt(4,x,d2);
   mpplt(5,x,hh);
   drawnow;
   [~,imx]=max(abs(d2));
   cp(k)=x(imx);
end
plot_fpmap(cp,flst,6);
fit_fpmap(cp,flst,pa);
for k=1:nfig
    print('-depsc2',sprintf('-f%d',k),sprintf('%s_%d.eps',fn,k));
end
return
end

%==============================================================

function mpfrm(n,lim1,lim2,lab)
mglab='magnitude (dB)';
phlab='phase (cyc)';
x_lab='distance from stapes (cm)';
figure(n);clf
subplot(211);hold on;axis(lim1);
title(lab);
ylabel(mglab);
subplot(212);hold on;axis(lim2);
xlabel(x_lab);
ylabel(phlab);
return
end

function mpplt(n,x,y)
figure(n);
ii=~isnan(y);
mg=20*log10(max(1e-9,abs(y)));
ph=unwrap(angle(y))/(2*pi);
subplot(211);plot(x(ii),mg(ii));
subplot(212);plot(x(ii),ph(ii));
return
end

function plot_fpmap(cp1,flst,fig)
fk=flst/1000;
cp0 = f2p(fk,0.1654,0.021,0.990,3.5);
figure(fig);clf
semilogx(fk,cp0,':k',fk,cp1,'-o');
xlabel('frequency (kHz)')
ylabel('distance from stapes (cm)')
axis([0.2 20 0 3])
fprintf('   f      x\n')
for k=1:length(flst)
    fprintf('%6.3f %6.3f\n',fk(k),cp1(k))
end
return
end

function p=f2p(fk,a,b,c,xl)
p = (1 - (log10((fk / a) + c) / b) / 100) * xl;
return
end

function fit_fpmap(cp1,flst,pa)
fk=flst/1000;
cp0=f2p(fk,0.1654,0.021,0.990,3.5); % man
c0=polyfit(log2(fk),cp0,1);
c1=polyfit(log2(fk),cp1,1);
pa=fitmap(c0,c1,pa);
wrtpar(pa,'fitmap.txt');
return
end

function pa=fitmap(c0,c1,pa)
dex=(c0(1)-c1(1))/2;
scl=2.^(c0(2)-c1(2));
fprintf('fitmap: dex=%6.3f scl=%6.3f\n',dex,scl);
pa.k1o=pa.k1o*scl;
pa.k2o=pa.k2o*scl;
pa.k3o=pa.k3o*scl;
pa.k4o=pa.k4o*scl;
pa.k5o=pa.k5o*scl;
pa.m1o=pa.m1o/scl;
pa.m2o=pa.m2o/scl;
pa.k1e=pa.k1e-dex;
pa.k2e=pa.k2e-dex;
pa.k3e=pa.k3e-dex;
pa.k4e=pa.k4e-dex;
pa.k5e=pa.k5e-dex;
pa.m1e=pa.m1e+dex;
pa.m2e=pa.m2e+dex;
return
end

%================================================

function pa=modpar22
pa.tl = 1;                  % transmission-line model
% default parameters
pa.aco = 0.01;              % area of cochlea
pa.ace = -0.4;              % cochlear-area taper
pa.ast = 0.01;              % area of stapes
pa.arw = 0.01;              % area of round window
pa.krw = 5e7;               % stiffness of round window
pa.khe = 0;                 % stiffness of helicotrema
pa.rhe = 0;                 % damping of helicotrema
pa.mhe = 1e-5;              % mass of helicotrema
pa.xhe = 0;                 % length of helicotrema
pa.ahe= 0;                  % area of helicotrema
pa.mmeq = 9;                % micromechanics eqs.
pa.dof = 2;                 % degrees-of-freedom
pa.gmo=1;                   % gamma at base [0.5]
pa.gme=0;                   % basal-gamma taper [-8]
pa.t1o=1e-9;                % OHC time constant
pa.t1e=0;                   % time-constant taper
pa.r1c=0;                   % 1st-DOF damping constant
pa.r2c=0;                   % 2nd-DOF damping constant
pa.m=1;                     % transmission-line model
pa.k5o= 1e5;
pa.k5e= 0;
pa.r4o= 0;
pa.r4e= 0;
pa.r5o= 400;
pa.r5e= 0;
%    fdm16 : model parameters for JASA manuscript
pa.n=3501;
pa.xl = 3.500;
pa.yw = 0.100;
pa.zh = 0.100;
pa.rho= 1.000;
pa.bwo= 0.050;
pa.bwe= 0.000;
pa.gpo= 1.000;
pa.gam= 1.000;
% partition impedance
pa.k1o=2.394e+08;
pa.k2o=3.036e+08;
pa.k3o=3.151e+08;
pa.k4o=4.045e+08;
pa.r1o=8.714e+02;
pa.r2o=1.979e+03;
pa.r3o=1.049e+00;
pa.m1o=7.417e-03;
pa.m2o=3.417e-02;
pa.k1e=-2.6655;
pa.k2e=-3.4762;
pa.k3e=-2.9092;
pa.k4e=-2.7948;
pa.r1e=-1.3937;
pa.r2e=-1.2466;
pa.r3e=0.1033;
pa.m1e=-0.0628;
pa.m2e=-0.0828;
pa.xtap=0.2239;        % taper at apex
return
end

%==============================================================

function [xx,Yb,Pd,Db,Dh]=fdmod22(pa,flst) % two chamber
N=pa.n;
nf=length(flst);
Yb=zeros(N,nf);
Pd=zeros(N,nf);
Db=zeros(N,nf);
Dh=zeros(N,nf);
% loop over frequencies
for k=1:nf
    f=flst(k);
    [a,q,x,y,D]=mxfill(f,pa);   % initialize arrays
    p=mxsolve(a,q);
    [pd,yb,dd]=p_dif(p,y,f,D);
    [pd,yb,dd]=p_nan(pd,yb,dd); % restrict plotting range
    Pd(:,k)=pd;
    Yb(:,k)=yb;
    Db(:,k)=dd(:,1);
    Dh(:,k)=dd(:,2);
end
xx=x/pa.xl;
return
end

%==============================================================

function [a,q,x,y,D]=mxfill(f,pa)
% set parameter values
chsz=pa.chsz;
chsz=chsz*2/sum(chsz); % normalize channel sizes;
m=length(chsz);
n=pa.n;
xl=pa.xl; rho=pa.rho; ast=pa.ast; ahe=pa.ahe; g=pa.gam;
aco=pa.aco; ace=pa.ace; bwo=pa.bwo; bwe=pa.bwe;
% useful constructs
s=2i*pi*f;
dx=xl/(n-1);
srd=s*rho*dx;
% compute admittance for all x
x=transpose(linspace(0,xl,n));
[z1,z2,zh,za]=imped(x,s,pa);
ac=aco*exp(ace*x);
bw=bwo*exp(bwe*x);
y=zeros(n,2,m);
A=zeros(n,m);
Y=zeros(n,m,m);
a=zeros(3,m,m,n);
% initialize partition admittance
D=[1 -1;0 0];
B=[-1;1];
H=[-1 -1;-1 1];
zg=zh-g*za;
for k=1:n
    A(k,:)=ac(k)*chsz;
    hh=z2(k)/(z2(k)+zh(k));
    zk=z1(k)+hh*zg(k);
    % expand zk
    zk(:,m)=0;
    zk(m,:)=1;
    Y(k,:,:)=zk\D;
    y(k,1,1)=1/zk(1);
    y(k,2,1)=hh/zk(1);
end
% initialize block-tridiagonal matrix a
for k=2:(n-1)
    A0=diag(A(k,:));
    A1=diag(A(k+1,:));
    yr=bw(k)*dx*squeeze(Y(k,:,:));
    a(1,:,:,k)=-A0;
    a(2,:,:,k)=A0+A1+srd*yr;
    a(3,:,:,k)=-A1;
end
% boundary conditions
A1=diag(A(1,:));
An=diag(A(n,:));
a(1,:,:,1)=zeros(m,m);
a(2,:,:,1)=2*A1;
a(3,:,:,1)=-2*A1;
a(1,:,:,n)=-An;
a(2,:,:,n)=An+ahe*H;
a(3,:,:,n)=zeros(m,m);
% initialize q
q=zeros(m,n);
q(:,1)=-2*s*srd*ast*B;
return
end

function [z1,z2,z3,z4]=imped(x,s,pa)
% compute impedance for all x
x=x.*(1+(pa.xtap*x).^6);
z1=pa.k1o*exp(pa.k1e*x)/s+pa.r1o*exp(pa.r1e*x)+pa.m1o*exp(pa.m1e*x)*s+pa.r1c;
z2=pa.k2o*exp(pa.k2e*x)/s+pa.r2o*exp(pa.r2e*x)+pa.m2o*exp(pa.m2e*x)*s+pa.r2c;
z3=pa.k3o*exp(pa.k3e*x)/s+pa.r3o*exp(pa.r3e*x);
z4=pa.k4o*exp(pa.k4e*x)/s+pa.r4o*exp(pa.r4e*x);
return
end

function p=mxsolve(a,q)
% solve tri-diagnonal block matrix equation
n=length(q);
% work down from top
b=squeeze(a(2,:,:,1));
a(3,:,:,1)=b\squeeze(a(3,:,:,1));
q(:,1)=b\q(:,1);
for k=2:n
   c=squeeze(a(1,:,:,k))*squeeze(a(3,:,:,k-1));
   p=squeeze(a(1,:,:,k))*q(:,k-1);
   b=squeeze(a(2,:,:,k))-c;
   a(3,:,:,k)=b\squeeze(a(3,:,:,k));
   q(:,k)=b\(q(:,k)-p);   
end
% work up from bottom
for k=(n-1):-1:1
   q(:,k)=q(:,k)-squeeze(a(3,:,:,k))*q(:,k+1);
end
p=q;
return
end

function [pd,yd,dd]=p_dif(p,y,f,D)
pref=2.848e-4; % dB SPL pressure reference (dyn/cm^2)
dref=1e-7;     % 1-nm displacement reference (cm)
n=size(p,2);
pd=zeros(n,1);
yd=zeros(n,1);
vd=zeros(n,2);
pdref=D(1,:)*p(:,1);
for k=1:n
   pk=D*p(:,k)/pdref;
   yk=squeeze(y(k,:,:));
   vk=yk*pk;
   vd(k,:)=vk(:)*(pref/dref);  % ST,HB velocity
   pd(k)=pk(1)-pk(end);        % ST-SV pressure difference
   yd(k)=vk(1)/pd(k);          % BM admittance
end
dd=vd./(2i*pi*f);              % BM,HB displacement
return
end

function [pd,yd,dd]=p_nan(pd,yd,dd)
% use nan to restrict plotted range
n=length(pd);
for k=1:n
   if(abs(dd(k,1))<1e-9)   
      ii=k:n;
      pd(ii,:)=nan;
      yd(ii,:)=nan;
      dd(ii,:)=nan;
      break;
   end
end   
return
end
