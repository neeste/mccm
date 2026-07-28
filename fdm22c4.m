% fdm22c4.m - compute frequency-domain four-chamber cochlear model 
function fdm22c4
ntc=1;
%
% compute 4ch cochlear model 
pa=modpar22;
% chamber sizes
pa.chsz=[1.1 0.3 0.1 0.5]; % err=28.14 28.14
pa.chpd=[1.000 -1.489 -0.589 -1.107 ; -0.408 1.000 0.251 -0.809 ; -0.611 -0.587 1.000 -0.859]; % err=28.08 28.02
pa.hbt=-20.7878;
if (ntc)
    flst=500*2.^(1:3);
    Dt=10.^(-tc_x(flst,pa.n,0.97,0)/20);
    Pt=[];
else
    [rt,ph,flst]=excpat79(pa.N);
    Dt=10.^((rt-30)/20);
    Pt=exp(2i*pi*ph);
end
% compute frequency-domain cochlear model
[xx,Yb,Pd,Db,Dh]=fdmod22(pa,flst);
% plot model results
fdm_plt(flst,Dt,Pt,xx,Yb,Pd,Db,Dh,pa.hbt);
end

function fdm_plt(flst,Dt,Pt,xx,Yb,Pd,Db,Dh,hbt)
Zb=1./Yb;
Hh=Dh./Db;
mpxplt(1,xx,Yb,[-100    0],[-0.5  0.5],'admittance');
rixplt(2,xx,Zb,[-250  250],[-250  250],'impedance');
mpxplt(3,xx,Pd,[ -60   40],[-7    1  ],'pressure re stapes');
mpxplt(4,xx,Db,[-100    0],[-7    1  ],'BM displacement (nm at 0 dB SPL)');
mpxplt(5,xx,Dh,[-120  -20],[-7    1  ],'HB displacement (nm at 0 dB SPL)');
mpxplt(6,xx,Hh,[ -50   40],[-0.5  1  ],'HB filter (dB)');
excpat(7,xx,Dh,[-100    0],[-7    1  ],'excitation pattern',hbt,Dt,Pt,flst);
drawnow
end

%================================================

function mpxplt(fig,xx,yy,dblim,phlim,lab)
mg=20*log10(abs(yy));
ph=angle(yy);
ph(ph<0)=ph(ph<0)+2*pi;
nf=size(yy,2);
for k=1:nf
    phk=unwrap(ph(:,k))/(2*pi);
    phr=max(phk)-min(phk);
    if (phr<1.01)
        phk=phk-round(mean(phk));
    end
    ph(:,k)=phk;
end
figure(fig);clf;
subplot(211)
plot(xx,mg)
axis([0 1 dblim])
title(lab)
ylabel('magnitude (dB)');
subplot(212)
plot(xx,ph)
axis([0 1 phlim]);
xlabel('distance from stapes (x/L)');
ylabel('phase (cyc)');
return
end

function rixplt(fig,xx,yy,rlim,ilim,lab)
z=zeros(size(xx));
rp=real(yy);
ip=imag(yy);
figure(fig);clf;
subplot(211)
set(gca,'ColorOrderIndex',1);
plot(xx,z,':k',xx,rp)
axis([0 1 rlim])
title(lab)
ylabel('real');
subplot(212)
set(gca,'ColorOrderIndex',1);
plot(xx,z,':k',xx,ip)
axis([0 1 ilim]);
xlabel('distance from stapes (x/L)');
ylabel('imaginary');
return
end

function excpat(fig,xx,Dh,dblim,phlim,lab,hbt,Dt,Pt,flst)
tpi=2*pi;
nf=length(flst);
db1=20*log10(abs(Dh))-hbt;
db2=20*log10(abs(Dt));
ph1=unwrap(angle(Dh))/tpi;
ph1(db1>120)=nan;
for k=1:nf
    if (ph1(1,k)<0)
        ph1(:,k)=ph1(:,k)+1;
    end
end
if (isempty(Pt))
    ph2=nan(size(Dt));
    for k=1:nf
        phd=ph1(:,k)-ph2(:,k);
        if (mean(phd(~isnan(phd)))>0.5)
            ph2(:,k)=ph2(:,k)+1;
        end
    end
    if (nf>2)
        if ((mean(ph1(1,:))>0)&&(ph1(1,1)<0))
            ph1(:,1)=ph1(:,1)+1;
        elseif ((mean(ph1(1,:))<0)&&(ph1(1,nf)>0))
            ph1(:,nf)=ph1(:,nf)-1;
        end
    end
else
    ph2=unwrap(angle(Pt))/tpi;
    n2=round(length(Dh)/2);
    n3=round(length(Dh)/3);
    for k=1:nf
        db2(:,k)=db2(:,k)+(db1(n2,k)-db2(n2,k));
        ph2(:,k)=ph2(:,k)+(ph1(n3,k)-ph2(n3,k));
    end
    phlim(1)=-5;
end
figure(fig);clf;
subplot(211);
plot(xx,db1); hold on
ax=gca;ax.ColorOrderIndex=1;
plot(xx,db2); hold off
axis([0 1 dblim]);
ylabel('magnitude (dB)');
title(lab);
subplot(212);plot(xx,ph1); hold on
set(gca,'ColorOrderIndex',1);
subplot(212);plot(xx,ph2); hold off
axis([0 1 phlim])
xlabel('x/L');
ylabel('phase (cyc)');
end

%==============================================================

function [xx,Yb,Pd,Db,Dh]=fdmod22(pa,flst)
nx=pa.n;
nf=length(flst);
Pd=zeros(nx,nf);
Yb=zeros(nx,nf);
Db=zeros(nx,nf);
Dh=zeros(nx,nf);
Hh=zeros(nx,nf);
for k=1:nf
   f=flst(k);
   [x,pd,yb,d1,d2,hh]=model(f,pa);
   Pd(:,k)=pd;
   Yb(:,k)=yb;
   Db(:,k)=d1;
   Dh(:,k)=d2;
   Hh(:,k)=hh;
end
xx=x/pa.xl;
end

function [x,pd,yb,d1,d2,hh]=model(f,pa)
% initialize arrays
[a,q,x,y,D]=mxfill(f,pa);
p=mxsolve(a,q);
[pd,yb,dd,hh]=p_dif(p,y,f,D);
[pd,yb,dd,hh]=p_nan(pd,yb,dd,hh);      % restrict plotting range
d1=dd(:,1); % BM displacement
d2=dd(:,2); % HB displacement
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
[z1,z2,zo,za]=imped(x,s,pa);
ac=aco*exp(ace*x);
bw=bwo*exp(bwe*x);
y=zeros(n,2,m);
A=zeros(n,m);
Y=zeros(n,m,m);
a=zeros(3,m,m,n);
% initialize partition admittance
D=[pa.chpd; 0 0 0 0];
B=[1; -1; 0; 0];
H=[1 0 0 -1; 0 0 0 0; 0 0 0 0; -1 0 0 1];
for k=1:n
    A(k,:)=ac(k)*chsz;
    zk=[z1(k)  0      0;
         0    z2(k)   0
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
end

function [z1,z2,z3,z4]=imped(x,s,pa)
% compute impedance for all x
q=x.^2;
k1=pa.k1o*exp(pa.k1e*x+pa.k1q*q);
r1=pa.r1o*exp(pa.r1e*x+pa.r1q*q);
m1=pa.m1o*exp(pa.m1e*x+pa.m1q*q);
k2=pa.k2o*exp(pa.k2e*x+pa.k2q*q);
r2=pa.r2o*exp(pa.r2e*x+pa.r2q*q);
m2=pa.m2o*exp(pa.m2e*x+pa.m2q*q);
k3=pa.k3o*exp(pa.k3e*x+pa.k3q*q);
r3=pa.r3o*exp(pa.r3e*x+pa.r3q*q);
m3=pa.m3o*exp(pa.m3e*x+pa.m3q*q);
k4=pa.k4o*exp(pa.k4e*x+pa.k4q*q);
r4=pa.r4o*exp(pa.r4e*x+pa.r4q*q);
z1=k1/s+r1+m1*s;
z2=k2/s+r2+m2*s;
z3=k3/s+r3+m3*s;
z4=k4/s+r4;
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
end

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
end

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
end

%================================================

function pa=modpar22
% misc constants
pa.g=1;
pa.b=0.4;
pa.gam=1;
pa.rho=1;
pa.hbt=-20;
%pa.N=1000;
% Middle-ear stifness, damping, and mass
pa.km=2.1e5;
pa.cm=400;
pa.mm=0.045;
pa.As=0.01;
pa.Am=0.25;
pa.Gm=0.5;
pa.Pe=2.848e-4;
% four-chamber cochlear model - 2022
pa.n=1001;
pa.xl = 2.500;
pa.yw = 0.100;
pa.zh = 0.100;
pa.rho= 1.000;
pa.bwo= 0.050;
pa.bwe= 0.000;
pa.gpo= 1.000;
pa.gam= 1.000;
pa.aco = 0.01;              % area of cochlea
pa.ace = -0.4;              % cochlear-area taper
pa.ast = 0.01;              % area of stapes
pa.arw = 0.01;              % area of round window
pa.ahe= 0;                  % area of helicotrema
% chamber sizes
pa.chsz=[1.1 0.3 0.1 0.5]; % err=28.14 28.14
pa.chpd=[1.000 -1.489 -0.589 -1.107 ; -0.408 1.000 0.251 -0.809 ; -0.611 -0.587 1.000 -0.859]; % err=28.08 28.02
% ------ parfit ------
pa.k1o=1.804e+08;
pa.r1o=108.6;
pa.m1o=5.653e-08;
pa.k2o=9.998e+09;
pa.r2o=8.936;
pa.m2o=0.002121;
pa.k3o=7.663e+04;
pa.r3o=0.001989;
pa.m3o=57.74;
pa.k4o=8.449e-06;
pa.r4o=1.401e-05;
pa.aco=0.003395;
pa.k1e=-3.1356;
pa.r1e=-0.8374;
pa.m1e=-0.1881;
pa.k2e=-6.5;
pa.r2e=0.1843;
pa.m2e=0.8440;
pa.k3e=0.0392;
pa.r3e=-0.0066;
pa.m3e=-0.0024;
pa.k4e=0.0061;
pa.r4e=-0.0193;
pa.ace=0.2890;
pa.k1q=0.000155;
pa.r1q=0.000213;
pa.m1q=-0.000016;
pa.k2q=0.7;
pa.r2q=0.000016;
pa.m2q=-0.000381;
pa.k3q=0.000013;
pa.r3q=0.000000;
pa.m3q=-0.000002;
pa.k4q=0.000134;
pa.r4q=0.000957;
pa.acq=0.000000;
pa.hbt=-20.7878;
% ------
end

