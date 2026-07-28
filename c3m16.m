% three-chamber cochlear model

function c3m16
fn='c3m16';nfig=0;
% fetch model parameters
pa=modpar16;
pa.m=3;
% prepare figures
L=pa.xl;
mpfrm(1,[0 L -100   0],[0 L -0.5 0.5],'BM admittance');
mpfrm(2,[0 L  -80  20],[0 L -9.0 1.0],'BM presure difference');
mpfrm(3,[0 L -110 -10],[0 L -9.0 1.0],'BM displacement re stapes');
mpfrm(4,[0 L -110 -10],[0 L -9.0 1.0],'HB displacement re stapes');
mpfrm(5,[0 L  -15   5],[0 L -0.1 0.3],'HB transfer');
% loop over 10 frequencies
flst=400*(sqrt(2).^(0:9));
nf=length(flst);
cp=zeros(1,nf);
for k=1:nf
   f=flst(k);
   [x,pd,yd,dd,hh]=model(f,pa);
   % plot data
   mpplt(1,x,yd(:,1));
   mpplt(2,x,pd(:,1));
   mpplt(3,x,dd(:,1));
   mpplt(4,x,dd(:,2));
   mpplt(5,x,hh(:));
   drawnow;
   [~,imx]=max(abs(dd(:,2)));
   cp(k)=x(imx);
end
plot_fpmap(cp,flst,6);
for k=1:nfig
    print('-depsc2',sprintf('-f%d',k),sprintf('%s_%d.eps',fn,k));
end
return

function plot_fpmap(cp,flst,fig)
fk=flst/1000;
a=0.1654; b=0.021; c=0.990; xl=3.5;
cp0 = (1 - (log10((fk / a) + c) / b) / 100) * xl;
figure(fig);clf
semilogx(fk,cp0,':k',fk,cp,'-o');
xlabel('frequency (kHz)')
ylabel('distance from stapes (cm)')
axis([0.2 20 0 3])
fprintf('   f      x\n')
for k=1:length(flst)
    fprintf('%6.3f %6.3f\n',fk(k),cp(k))
end
return

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

function mpplt(n,x,y)
figure(n);
ii=~isnan(y);
mg=20*log10(max(1e-9,abs(y)));
ph=unwrap(angle(y))/(2*pi);
subplot(211);plot(x(ii),mg(ii));
subplot(212);plot(x(ii),ph(ii));
return

%==============================================================

function [x,pd,yd,dd,hh]=model(f,pa)
% initialize arrays
[a,q,x,y,D]=mxfill(f,pa);
p=mxsolve(a,q);
[pd,yd,dd,hh]=p_dif(p,y,f,D);
return

%==============================================================

function [a,q,x,y,D]=mxfill(f,pa)
% set parameter values
m=pa.m; n=pa.n;
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
zg=zh-g*za;
y=zeros(n,2,m);
A=zeros(n,m);
Y=zeros(n,m,m);
a=zeros(3,m,m,n);
% initialize partition admittance
chsz=[1 1 1e9];
D=[1 0 -1;0 0 0;0 0 0];
B=[-1;0;1];
H=[1 0 -1;-1 0 1;0 0 0];
for k=1:n
    zk=[z1(k)+zg(k) -zg(k);
        -zh(k) z2(k)+zh(k)];
    yp=zk\[1;0];
    zk(:,m)=0;
    zk(m,:)=1;
    A(k,:)=ac(k)*chsz;
    Y(k,:,:)=zk\D;
    y(k,1,1)=yp(1);
    y(k,2,1)=yp(1)-yp(2);
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

function [z1,z2,z3,z4]=imped(x,s,pa)
% compute impedance for all x
x=x.*(1+(pa.xtap*x).^6);
z1=pa.k1o*exp(pa.k1e*x)/s+pa.r1o*exp(pa.r1e*x)+pa.m1o*exp(pa.m1e*x)*s+pa.r1c;
z2=pa.k2o*exp(pa.k2e*x)/s+pa.r2o*exp(pa.r2e*x)+pa.m2o*exp(pa.m2e*x)*s+pa.r2c;
z3=pa.k3o*exp(pa.k3e*x)/s+pa.r3o*exp(pa.r3e*x);
z4=pa.k4o*exp(pa.k4e*x)/s+pa.r4o*exp(pa.r4e*x);
return

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

function [pd,yd,dd,hh]=p_dif(p,y,f,D)
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
   pd(k,:)=pk(1);
   yd(k,1)=vk(1)/pk(1);
   vd(k,:)=vk(:)*(pref/dref);
end
hh=vd(:,2)./vd(:,1);
dd=vd./(2i*pi*f);
% restrict plotted range
for k=1:n
   if(abs(vd(k,1))<1e-8)   
      ii=k:n;
      pd(ii,:)=nan;
      dd(ii,:)=nan;
      hh(ii,:)=nan;
      break;
   end
end   
return
