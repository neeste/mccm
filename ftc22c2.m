% ftc22c2.m - frequency-domain two-chamber cochlear model tuning curves
function ftc22c2
%
% compute 3DOF cochlear model 
pa=modpar22;
% chamber sizes
pa.chsz=[1 1]; % relative chamber sizes
flst=500*2.^(-1:5);
Dt=10.^(-tc_x(flst,pa.n,pa.xp,0)/20);
Pt=[];
% compute frequency-domain cochlear model
[xx,~,~,~,Dh]=fdmod22(pa,flst);
% locate Dh peaks
bp=best_place(flst,xx*pa.xl,Dh);
% plot model results
hbt=pa.hbt;
excpat(1,xx,Dh,[-100    0],[-7    1  ],'excitation pattern',hbt,Dt,Pt,flst);
fpmplt(2,bp,flst,pa.xl);
% compute model frequency-tuning curves
[f,epl]=tc_f(flst); % EPL average of 6 cats
spl=fdmftc(pa,f,bp);
ftcplt(3,f,epl,spl);
end

%================================================

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

%==============================================================

function fpmplt(fig,bp,flst,L)
fk=flst/1000;
xp=0.97; % portion of x mapped
%cp = f2p(fk,0.1654,0.021,0.99,3.5); % man
cp = f2p(fk,0.4560,0.021,0.80,2.5*xp); % cat
figure(fig);clf
ax=gca;ax.ColorOrderIndex=1;
semilogx(fk,cp,'--k',fk,bp,'-o');
xlabel('frequency (kHz)')
ylabel('distance from stapes (cm)')
title('frequency-place map (cat)')
axis([0.2 20 0 L])
drawnow
end

function p=f2p(fk,a,b,c,xl)
p = (1 - (log10((fk / a) + c) / b) / 100) * xl;
return
end

%================================================

function ftcplt(fig,f,ftc,spl)
figure(fig);clf
ax=gca;ax.ColorOrderIndex=1;
semilogx(f,ftc,'--'); hold on
ax=gca;ax.ColorOrderIndex=1;
semilogx(f,spl); hold off
axis([0.1 20 0 90])
title('neural threshold frequency-tuning curves (EPL,cat)')
xlabel('frequency (kHz)')
ylabel('SPL (dB)')
text(10,5,'ftc22c2')
drawnow
end

%==============================================================

function spl=fdmftc(pa,flst,xlst)
flst=flst*1000;
spl=zeros(length(flst),length(xlst));
[xx,~,~,~,Dh]=fdmod22(pa,flst);
for k=1:length(xlst)
    [~,j]=min(abs(xlst(k)-xx(:)*pa.xl));
    spl(:,k)=pa.hbt-20*log10(abs(Dh(j,:)))';
end
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
pa.hbt= 20; pa.xp=0.97;
pa.xtap=0.0; pa.xtex=6;
pa.ast= 0.01;   % area of stapes
pa.khe=1.525; pa.rhe=1.096e-09; pa.mhe=55.89928; % err=12.73 12.73
pa.kme=8.066e+06; pa.rme=567.18; pa.mme=0.108456; % err=12.73 12.71
% ---- parfit ----
pa.k1o=3.295e+05;
pa.r1o=0.02935;
pa.m1o=1.355e-09;
pa.k2o=3.376e+06;
pa.r2o=9.314;
pa.m2o=3.913e-05;
pa.k3o=5.157e+05;
pa.r3o=5.255e-05;
pa.k4o=6.104e+05;
pa.r4o=1e-09;
pa.aco=0.01;
pa.k1e=-2.4151;
pa.r1e=3.2511;
pa.m1e=0.0013;
pa.k2e=-6.6464;
pa.r2e=-5.2152;
pa.m2e=-2.0247;
pa.k3e=-3.8631;
pa.r3e=2.0412;
pa.k4e=-3.7895;
pa.r4e=0.0000;
pa.ace=6.8386;
pa.k1q=0.187921;
pa.r1q=-1.255692;
pa.m1q=-0.000173;
pa.k2q=-0.000000;
pa.r2q=1.530265;
pa.m2q=0.597474;
pa.k3q=-0.354283;
pa.r3q=0.590922;
pa.k4q=0.763145;
pa.r4q=0.000000;
pa.acq=-1.922025;
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
Zh=pa.khe/s+pa.rhe+pa.mhe*s;
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
zz(n)=Zh;
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
