% fdm22.m - frequency-domain two-chamber cochlear model 
function fdm22c2
ntc=1;
%
% compute 3DOF cochlear model 
pa=modpar22;
% chamber sizes
pa.chsz=[1 1]; % relative chamber sizes
if (ntc)
    flst=500*2.^(-1:6);
    Dt=10.^(-tc_x(flst,pa.n)/20);
    Pt=[];
else
    [rt,ph,flst]=excpat79(pa.n);
    Dt=10.^((rt-30)/20);
    Pt=exp(2i*pi*ph);
end
% compute frequency-domain cochlear model
[xx,Yb,Pd,Db,Dh]=fdmod22(pa,flst);
% plot model results
fdm_plt(flst,Dt,Pt,xx,Yb,Pd,Db,Dh,pa.hbt);
end

%================================================

function fdm_plt(flst,Dt,Pt,xx,Yp,Pd,Db,Dh,hbt)
Zp=1./Yp;
Hh=Dh./Db;
mpxplt(1,xx,Yp,[ -80   20],[ -0.5  0.5],'CP admittance');
rixplt(2,xx,Zp,[  -5    5],[ -5    5  ],'CP impedance');
mpxplt(3,xx,Pd,[ -60   40],[ -7    1  ],'CP pressure difference');
mpxplt(4,xx,Db,[ -80   20],[ -7    1  ],'BM displacement (nm at 0 dB SPL)');
mpxplt(5,xx,Dh,[ -80   20],[ -7    1  ],'HB displacement (nm at 0 dB SPL)');
mpxplt(6,xx,Hh,[ -50   30],[ -0.5  0.5],'HB/BM filter (dB)');
excpat(7,xx,Dh,[-100    0],[ -7    1  ],'excitation pattern',hbt,Dt,Pt,flst);
drawnow
end

function mpxplt(fig,xx,yy,dblim,phlim,lab)
mg=20*log10(abs(yy));
ph=angle(yy);
ph(ph<0)=ph(ph<0)+2*pi;
nf=size(yy,2);
for k=1:nf
    phk=unwrap(ph(:,k))/(2*pi);
    phr=max(phk)-min(phk);
    if (phr<1.01)
        phk=phk-(max(phk)+min(phk))/2;
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
ip(isnan(rp))=nan;
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
if (isempty(Pt))
    ph2=nan(size(Dt));
else
    ph2=unwrap(angle(Pt))/tpi;
    for k=1:nf
        phd=ph1(:,k)-ph2(:,k);
        if (mean(phd(~isnan(phd)))>0.5)
            ph2(:,k)=ph2(:,k)+1;
        end
    end
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
ax=gca;ax.ColorOrderIndex=1;
subplot(212);plot(xx,ph2); hold off
axis([0 1 phlim])
xlabel('x/L');
ylabel('phase (cyc)');
end

%================================================

function pa=modpar22
% default parameters
pa.gam= 1.000;
pa.n=2501;
pa.xl = 2.500;
pa.yw = 0.100;
pa.zh = 0.100;
pa.rho= 1.000;
pa.bwo= 0.050;
pa.bwe= 0.000;
pa.hbt= 20;
pa.xtap=0.0; pa.xtex=6;
pa.ast= 0.01;   % area of stapes
pa.rhe=1.065e-09; pa.mhe=20.164338; % err=16.35 16.35
pa.kme=6.288e+06; pa.rme=528.94; pa.mme=0.060389; % err=16.35 16.35
% ---- parfit ----
pa.k1o=5.354e+05;
pa.r1o=0.06762;
pa.m1o=2.394e-09;
pa.k2o=2.526e+06;
pa.r2o=7.727;
pa.m2o=2.835e-05;
pa.k3o=6.166e+05;
pa.r3o=3.604e-05;
pa.k4o=7.694e+05;
pa.r4o=1e-09;
pa.aco=0.01;
pa.k1e=-2.3789;
pa.r1e=2.8601;
pa.m1e=0.0016;
pa.k2e=-6.8895;
pa.r2e=-5.4620;
pa.m2e=-1.6445;
pa.k3e=-3.7605;
pa.r3e=1.9455;
pa.k4e=-3.8173;
pa.r4e=0.0000;
pa.ace=5.6696;
pa.k1q=0.100895;
pa.r1q=-1.149295;
pa.m1q=-0.000759;
pa.k2q=-0.000006;
pa.r2q=1.653353;
pa.m2q=0.506121;
pa.k3q=-0.454258;
pa.r3q=0.741516;
pa.k4q=0.775636;
pa.r4q=0.000000;
pa.acq=-1.458174;
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

function [a,q,x,y,D]=mxfill(f,pa)
% set parameter values
chsz=pa.chsz;
chsz=chsz*2/sum(chsz); % normalize channel sizes;
m=length(chsz);
n=pa.n;
xl=pa.xl; rho=pa.rho; ast=pa.ast; g=pa.gam;
bwo=pa.bwo; bwe=pa.bwe;
% useful constructs
s=2i*pi*f;
dx=xl/(n-1);
srd=s*rho*dx;
Zm=pa.kme/s+pa.rme+pa.mme*s;
% compute admittance for all x
x=transpose(linspace(0,xl,n));
[z1,z2,zh,za,ac]=imped(x,s,pa);
bw=bwo*exp(bwe*x);
y=zeros(n,2,m);
A=zeros(n,m);
Y=zeros(n,m,m);
a=zeros(3,m,m,n);
% initialize partition admittance
D=[1 -1;0 0];
B=[-1;1];
zg=zh-g*za;
hh=z2./(z2+zh);
zz=z1+hh.*zg;
%nhe=round(n*0.99);
zhe=pa.rhe+pa.mhe*s;
zz(n)=zhe;
for k=1:n
    A(k,:)=ac(k)*chsz;
    zk=zz(k);
    % expand zk
    zk(:,m)=0;
    zk(m,:)=1;
    Y(k,:,:)=zk\D;
    y(k,1,1)=1/z2(k);
    y(k,2,1)=1/zz(k);
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
a(2,:,:,n)=An;
a(3,:,:,n)=zeros(m,m);
% initialize q
q=zeros(m,n);
q(:,1)=-2*(srd*ast./Zm)*B;
return
end

function [z1,z2,z3,z4,ac]=imped(x,s,pa)
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
ac=pa.aco*exp(pa.ace*x+pa.acq*q);
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

function [pd,yp,dd]=p_dif(p,y,f,D)
pref=2.848e-4; % dB SPL pressure reference (dyn/cm^2)
dref=1e-7;     % 1-nm displacement reference (cm)
n=size(p,2);
pd=zeros(n,1);
yp=zeros(n,1);
vd=zeros(n,2);
%pdref=D(1,:)*p(:,1);
pdref=1;
for k=1:n
   pk=D*p(:,k)/pdref;
   yk=squeeze(y(k,:,:));
   vk=yk*pk;
   vd(k,:)=vk(:)*(pref/dref);  % ST,HB velocity
   pd(k)=pk(1)-pk(end);        % ST-SV pressure difference
   yp(k)=vk(2)/pd(k);          % CP admittance
end
dd=vd./(2i*pi*f);              % convert velecity to displacement
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
