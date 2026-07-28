% c4m22 four-chamber cochlear model - 2022
function c4m22
fn='c4m22';nfig=6;
% fetch model parameters
pa=modpar22;
pa.chsz=[1 1]; % relative chamber sizes
pa.chsz=[1.205610 0.340113 0.003639 0.450637]; % err=1.4371
pa.chpd=[1.000 -1.489 -0.588 -1.105 ; -0.403 1.000 0.249 -0.808 ; -0.610 -0.585 1.000 -0.862]; % err=1.4381
pa.k5o=3.388e+04; pa.k5e=-0.010; pa.r5o=1.982e+00; pa.r5e=0.006; % err=1.4382
%pa=parfit22(pa,1);
% compute model results for frequency set
flst=400*(sqrt(2).^(0:9));
[x,pdf,ybf,d1f,d2f,hhf]=modset(flst,pa);
% prepare figures
L=pa.xl;
mpfrm(1,[0 L -100   0],[0 L -0.3 0.3],'BM admittance');
mpfrm(2,[0 L  -80  20],[0 L -9.0 1.0],'BM presure difference');
mpfrm(3,[0 L -110 -10],[0 L -9.0 1.0],'BM displacement re stapes');
mpfrm(4,[0 L -110 -10],[0 L -9.0 1.0],'HB displacement re stapes');
mpfrm(5,[0 L  -15  15],[0 L -0.3 0.3],'HB transfer');
% loop over 10 frequencies
nf=length(flst);
cp=zeros(1,nf);
for k=1:nf
   pd=pdf(:,k);
   yb=ybf(:,k);
   d1=d1f(:,k);
   d2=d2f(:,k);
   hh=hhf(:,k);
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

function p=f2p(fk,a,b,c,xl)
p = (1 - (log10((fk / a) + c) / b) / 100) * xl;
return

function fit_fpmap(cp1,flst,pa)
fk=flst/1000;
cp0=f2p(fk,0.1654,0.021,0.990,3.5); % man
c0=polyfit(log2(fk),cp0,1);
c1=polyfit(log2(fk),cp1,1);
pa=fitmap(c0,c1,pa);
wrtpar(pa,'fitmap.txt');
return

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

%==============================================================

function [x,pdf,ybf,d1f,d2f,hhf]=modset(flst,pa)
nx=pa.n;
nf=length(flst);
pdf=zeros(nx,nf);
ybf=zeros(nx,nf);
d1f=zeros(nx,nf);
d2f=zeros(nx,nf);
hhf=zeros(nx,nf);
for k=1:nf
   f=flst(k);
   [x,pd,yb,d1,d2,hh]=model(f,pa);
   pdf(:,k)=pd;
   ybf(:,k)=yb;
   d1f(:,k)=d1;
   d2f(:,k)=d2;
   hhf(:,k)=hh;
end
return

function [x,pd,yb,d1,d2,hh]=model(f,pa)
% initialize arrays
[a,q,x,y,D]=mxfill(f,pa);
p=mxsolve(a,q);
[pd,yb,dd,hh]=p_dif(p,y,f,D);
[pd,yb,dd,hh]=p_nan(pd,yb,dd,hh);      % restrict plotting range
d1=dd(:,1); % BM displacement
d2=dd(:,2); % HB displacement
return

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
[z1,z2,zh,za,zo]=imped(x,s,pa);
ac=aco*exp(ace*x);
bw=bwo*exp(bwe*x);
y=zeros(n,2,m);
A=zeros(n,m);
Y=zeros(n,m,m);
a=zeros(3,m,m,n);
% initialize partition admittance
if (m==2)
    D=[1 -1;0 0];
    B=[-1;1];
    H=[-1 -1;-1 1];
    hh=z2./(z2+zh);
    zs=(zh-g*za).*hh;
    for k=1:n
        A(k,:)=ac(k)*chsz;
        zk=z1(k)+zs(k);
        % expand zk
        zk(:,m)=0;
        zk(m,:)=1;
        % compute admittances
        Y(k,:,:)=zk\D;           % radial constraints
        y(k,1,:)=Y(k,1,:);       % BM admittance
        y(k,2,:)=hh(k)*Y(k,1,:); % HB admittance
    end
else
    D=[pa.chpd; 0 0 0 0];
    B=[1; -1; 0; 0];
    H=[1 0 0 -1; 0 0 0 0; 0 0 0 0; -1 0 0 1];
    zs=z2+zh;
    for k=1:n
        A(k,:)=ac(k)*chsz;
        zk=[z1(k)  0      0;
             0    zs(k)   0
             0     0    zo(k)];
        % expand zk
        zk(:,m)=0;
        zk(m,:)=1;
        % compute admittances
        Yr=zk\D;
        Za=[0 0 0 0; 0 0 0 0; 0 -g*za(k) 0 0; 0 0 0 0];
        Y(k,:,:)=(eye(m)-Yr*Za)\Yr;
        y(k,1,:)=Y(k,1,:);   % BM admittance
        y(k,2,:)=-Y(k,2,:);  % HB admittance
    end
end
%y=Y(:,1:2,:); % extract BM & HB admittances
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

function [z1,z2,z3,z4,z5]=imped(x,s,pa)
% compute impedance for all x
x=x.*(1+(pa.xtap*x).^6);
z1=pa.k1o*exp(pa.k1e*x)/s+pa.r1o*exp(pa.r1e*x)+pa.m1o*exp(pa.m1e*x)*s+pa.r1c;
z2=pa.k2o*exp(pa.k2e*x)/s+pa.r2o*exp(pa.r2e*x)+pa.m2o*exp(pa.m2e*x)*s+pa.r2c;
z3=pa.k3o*exp(pa.k3e*x)/s+pa.r3o*exp(pa.r3e*x);
z4=pa.k4o*exp(pa.k4e*x)/s+pa.r4o*exp(pa.r4e*x);
z5=pa.k5o*exp(pa.k5e*x)/s+pa.r5o*exp(pa.r5e*x);
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
   vd(k,:)=vk(:)*(pref/dref);  % ST,HB velocity
   pd(k)=pk(1)-pk(end);        % ST-SV pressure difference
   yd(k)=vk(1)/pd(k);          % BM admittance
end
hh=vd(:,2)./vd(:,1);           % HB/BM transfer
dd=vd./(2i*pi*f);              % BM,HB displacement
return

function [pd,yd,dd,hh]=p_nan(pd,yd,dd,hh)
% use nan to restrict plotted range
n=length(pd);
for k=1:n
   if(abs(dd(k,1))<1e-9)   
      ii=k:n;
      pd(ii,:)=nan;
      yd(ii,:)=nan;
      dd(ii,:)=nan;
      hh(ii,:)=nan;
      break;
   end
end   
return
