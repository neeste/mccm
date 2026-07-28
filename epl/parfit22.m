% parfit22.m - fit 3DOF model parameters to AN threshold
function parfit22
ntc=1;
%
% compute 3DOF cochlear model 
pa=modpar22;
if (ntc)
    flst=500*2.^(0:4);
    Dt=10.^(-tc_x(flst,pa.n)/20);
    Pt=[];
else
    [rt,ph,flst]=excpat79(pa.n);
    Dt=10.^((rt-30)/20);
    Pt=exp(2i*pi*ph);
end
for k=1:1000
    pa=fitpar22(pa,flst,Dt,Pt,0);
    fdm_plt(pa,flst,Dt,Pt);
end
end

function fdm_plt(pa,flst,Dt,Pt)
% compute frequency-domain cochlear model
[xx,Yb,Pd,Db,Dh,Hh]=modset(flst,pa);
% plot model results
Zb=1./Yb;
mpxplt(1,xx,Yb,[-100   0],[-0.5 0.5],'admittance');
rixplt(2,xx,Zb,[-100 200],[-200 100],'impedance');
mpxplt(3,xx,Pd,[ -60   40],   [-7 1],'pressure re stapes');
mpxplt(4,xx,Db,[-120  -20],   [-7 1],'BM displacement (nm at 0 dB SPL)');
mpxplt(5,xx,Dh,[-120  -20],   [-7 1],'HB displacement (nm at 0 dB SPL)');
mpxplt(6,xx,Hh,[ -50   30],   [-1 1],'HB filter (dB)');
excpat(7,xx,Dh,[-100    0],   [-7 1],'excitation pattern',pa,Dt,Pt,flst);
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
figure(fig);clf;
subplot(211)
plot(xx,z,':k',xx,rp)
axis([0 1 rlim])
title(lab)
ylabel('real');
subplot(212)
plot(xx,z,':k',xx,ip)
axis([0 1 ilim]);
xlabel('distance from stapes (x/L)');
ylabel('imaginary');
return
end

function excpat(fig,xx,Dh,dblim,phlim,lab,pa,Dt,Pt,flst)
tpi=2*pi;
nf=length(flst);
db1=20*log10(abs(Dh))-pa.hbt;
db2=20*log10(abs(Dt));
ph1=unwrap(angle(Dh))/tpi;
ph1(db1>120)=nan;
if (isempty(Pt))
    ph2=nan(size(Dt));
else
    ph2=unwrap(angle(Pt))/tpi;
end
for k=1:nf
    phd=ph1(:,k)-ph2(:,k);
    if (mean(phd(~isnan(phd)))>0.5)
        ph2(:,k)=ph2(:,k)+1;
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

function pa=fitpar22(pa,flst,Dt,Pt,ot)
if (ot==0)
    pv=getpar(pa);
elseif (ot==1)
    pv=pa.chsz(1:3);
elseif (ot==2)
    pv=pa.chpd(:);
end
% fit parameters to neural data
e1=refdev(pv,pa,flst,Dt,Pt,ot);
op=optimset('MaxFunEvals',1000);
pv=fminsearch(@(pv) refdev(pv,pa,flst,Dt,Pt,ot),pv,op);
e2=refdev(pv,pa,flst,Dt,Pt,ot);
pa=modpar22;
if (ot==0)
    pa=setpar(pa,pv);
    prnpar(pa);
elseif (ot==1)
    pa.chsz=[pv 2-sum(pv)];
    fprintf('pa.chsz=[%.6f %.6f %.6f %.6f]; %% ',pa.chsz)
elseif (ot==2)
    chpd=[pv(1:3),pv(4:6),pv(7:9),pv(10:12)];
    pa.chpd=[chpd(1,:)/chpd(1,1);chpd(2,:)/chpd(2,2);chpd(3,:)/chpd(3,3)];
    fprintf('pa.chpd=[%.3f %.3f %.3f %.3f ; ',pa.chpd(1,:));
    fprintf('%.3f %.3f %.3f %.3f ; ',         pa.chpd(2,:));
    fprintf('%.3f %.3f %.3f %.3f]; %% ',      pa.chpd(3,:));
end
fprintf('err=%.4g %.4g\n',e1,e2)
end

function err=refdev(pv,pa,flst,Dt,Pt,ot)
err=0;
if (ot==0)
    no=(length(pv)-1)/3; io=1:no; ie=io+no; % pv indices
    if (min(pv(io))<1e-9) err=1e3; return; end
    err=err+mean(abs(pv(ie)));
    r1=pv(2)*exp(pv(15));
    if (r1<90) err=err+90-r1; end
    pa=setpar(pa,pv);
elseif (ot==1)
    if (min(pv)<1e-3) err=1e3; return; end
    pa.chsz=[pv 2-sum(pv)];
elseif (ot==2)
    chpd=[pv(1:3),pv(4:6),pv(7:9),pv(10:12)];
    pa.chpd=[chpd(1,:)/chpd(1,1);chpd(2,:)/chpd(2,2);chpd(3,:)/chpd(3,3)];
end
% compute fdm22
[~,Yb,~,~,Dh]=modset(flst,pa);
D1=20*log10(abs(Dh))-pa.hbt;
D2=20*log10(abs(Dt));
DD=D1-D2;
[nx,nf]=size(DD);
if (flst(1)==500)
    ph1=unwrap(angle(Dh(:,:)))/(2*pi);
    phd=(ph1(1,:)-ph1(nx,:));
    err=err+abs(min(phd)-2.5)*4;
end
for k=1:nf
    pk=max(abs(diff(diff(log(abs(Yb(:,k)))))))*nx/2;
    err=err+pk;
    ii=~isnan(DD(:,k))&(D1(:,k)>-150);
    err=err+mean(abs(DD(ii,k)).^2);
    if (isempty(Pt))
        ph=unwrap(angle(Dh(:,k)))/(2*pi);
        ep=ph(end);
        tp=-3.6; % target phase
        if (ep>tp)
            err=err+(ep-tp)*10;
        end
    else
        PP=unwrap(angle(Dh(:,k)./Pt(:,k)));
        ii=~isnan(PP);
        err=err+mean(abs(PP(ii)).^2)*20;
    end
end
err=err/nf;
end

function pv=getpar(pa)
k=0;
k=k+1;pv(k)=pa.k1o;
k=k+1;pv(k)=pa.r1o;
k=k+1;pv(k)=pa.m1o;
k=k+1;pv(k)=pa.k2o;
k=k+1;pv(k)=pa.r2o;
k=k+1;pv(k)=pa.m2o;
k=k+1;pv(k)=pa.k3o;
k=k+1;pv(k)=pa.r3o;
k=k+1;pv(k)=pa.m3o;
k=k+1;pv(k)=pa.k4o;
k=k+1;pv(k)=pa.r4o;
k=k+1;pv(k)=pa.aco;
k=k+1;pv(k)=pa.k1e;
k=k+1;pv(k)=pa.r1e;
k=k+1;pv(k)=pa.m1e;
k=k+1;pv(k)=pa.k2e;
k=k+1;pv(k)=pa.r2e;
k=k+1;pv(k)=pa.m2e;
k=k+1;pv(k)=pa.k3e;
k=k+1;pv(k)=pa.r3e;
k=k+1;pv(k)=pa.m3e;
k=k+1;pv(k)=pa.k4e;
k=k+1;pv(k)=pa.r4e;
k=k+1;pv(k)=pa.ace;
k=k+1;pv(k)=pa.k1q;
k=k+1;pv(k)=pa.r1q;
k=k+1;pv(k)=pa.m1q;
k=k+1;pv(k)=pa.k2q;
k=k+1;pv(k)=pa.r2q;
k=k+1;pv(k)=pa.m2q;
k=k+1;pv(k)=pa.k3q;
k=k+1;pv(k)=pa.r3q;
k=k+1;pv(k)=pa.m3q;
k=k+1;pv(k)=pa.k4q;
k=k+1;pv(k)=pa.r4q;
k=k+1;pv(k)=pa.acq;
k=k+1;pv(k)=pa.hbt;
end

function pa=setpar(pa,pv)
k=0;
k=k+1;pa.k1o=pv(k);
k=k+1;pa.r1o=pv(k);
k=k+1;pa.m1o=pv(k);
k=k+1;pa.k2o=pv(k);
k=k+1;pa.r2o=pv(k);
k=k+1;pa.m2o=pv(k);
k=k+1;pa.k3o=pv(k);
k=k+1;pa.r3o=pv(k);
k=k+1;pa.m3o=pv(k);
k=k+1;pa.k4o=pv(k);
k=k+1;pa.r4o=pv(k);
k=k+1;pa.aco=pv(k);
k=k+1;pa.k1e=pv(k);
k=k+1;pa.r1e=pv(k);
k=k+1;pa.m1e=pv(k);
k=k+1;pa.k2e=pv(k);
k=k+1;pa.r2e=pv(k);
k=k+1;pa.m2e=pv(k);
k=k+1;pa.k3e=pv(k);
k=k+1;pa.r3e=pv(k);
k=k+1;pa.m3e=pv(k);
k=k+1;pa.k4e=pv(k);
k=k+1;pa.r4e=pv(k);
k=k+1;pa.ace=pv(k);
k=k+1;pa.k1q=pv(k);
k=k+1;pa.r1q=pv(k);
k=k+1;pa.m1q=pv(k);
k=k+1;pa.k2q=pv(k);
k=k+1;pa.r2q=pv(k);
k=k+1;pa.m2q=pv(k);
k=k+1;pa.k3q=pv(k);
k=k+1;pa.r3q=pv(k);
k=k+1;pa.m3q=pv(k);
k=k+1;pa.k4q=pv(k);
k=k+1;pa.r4q=pv(k);
k=k+1;pa.acq=pv(k);
k=k+1;pa.hbt=pv(k);
end

function prnpar(pa)
fp=fopen('modpar22.txt','wt');
fprintf(fp,'pa.k1o=%.4g;\n',pa.k1o);
fprintf(fp,'pa.r1o=%.4g;\n',pa.r1o);
fprintf(fp,'pa.m1o=%.4g;\n',pa.m1o);
fprintf(fp,'pa.k2o=%.4g;\n',pa.k2o);
fprintf(fp,'pa.r2o=%.4g;\n',pa.r2o);
fprintf(fp,'pa.m2o=%.4g;\n',pa.m2o);
fprintf(fp,'pa.k3o=%.4g;\n',pa.k3o);
fprintf(fp,'pa.r3o=%.4g;\n',pa.r3o);
fprintf(fp,'pa.m3o=%.4g;\n',pa.m3o);
fprintf(fp,'pa.k4o=%.4g;\n',pa.k4o);
fprintf(fp,'pa.r4o=%.4g;\n',pa.r4o);
fprintf(fp,'pa.aco=%.4g;\n',pa.aco);
fprintf(fp,'pa.k1e=%.4f;\n',pa.k1e);
fprintf(fp,'pa.r1e=%.4f;\n',pa.r1e);
fprintf(fp,'pa.m1e=%.4f;\n',pa.m1e);
fprintf(fp,'pa.k2e=%.4f;\n',pa.k2e);
fprintf(fp,'pa.r2e=%.4f;\n',pa.r2e);
fprintf(fp,'pa.m2e=%.4f;\n',pa.m2e);
fprintf(fp,'pa.k3e=%.4f;\n',pa.k3e);
fprintf(fp,'pa.r3e=%.4f;\n',pa.r3e);
fprintf(fp,'pa.m3e=%.4f;\n',pa.m3e);
fprintf(fp,'pa.k4e=%.4f;\n',pa.k4e);
fprintf(fp,'pa.r4e=%.4f;\n',pa.r4e);
fprintf(fp,'pa.ace=%.4f;\n',pa.ace);
fprintf(fp,'pa.k1q=%.6f;\n',round(pa.k1q,6));
fprintf(fp,'pa.r1q=%.6f;\n',round(pa.r1q,6));
fprintf(fp,'pa.m1q=%.6f;\n',round(pa.m1q,6));
fprintf(fp,'pa.k2q=%.6f;\n',round(pa.k2q,6));
fprintf(fp,'pa.r2q=%.6f;\n',round(pa.r2q,6));
fprintf(fp,'pa.m2q=%.6f;\n',round(pa.m2q,6));
fprintf(fp,'pa.k3q=%.6f;\n',round(pa.k3q,6));
fprintf(fp,'pa.r3q=%.6f;\n',round(pa.r3q,6));
fprintf(fp,'pa.m3q=%.6f;\n',round(pa.m3q,6));
fprintf(fp,'pa.k4q=%.6f;\n',round(pa.k4q,6));
fprintf(fp,'pa.r4q=%.6f;\n',round(pa.r4q,6));
fprintf(fp,'pa.acq=%.6f;\n',round(pa.acq,6));
fprintf(fp,'pa.hbt=%.4f;\n',pa.hbt);
fclose(fp);
end

%================================================

function pa=modpar22
% misc constants
pa.g=1;
pa.b=0.4;
pa.gam=1;
pa.rho=1;
% Middle-ear stifness, damping, and mass
pa.km=2.1e5;
pa.cm=400;
pa.mm=0.045;
pa.As=0.01;
pa.Am=0.25;
pa.Gm=0.5;
pa.Pe=2.848e-4;
% four-chamber cochlear model - 2022
pa.n=2501;
pa.xl = 2.500;
pa.yw = 0.100;
pa.zh = 0.100;
pa.rho= 1.000;
pa.bwo= 0.050;
pa.bwe= 0.000;
pa.gpo= 1.000;
pa.gam= 1.000;
pa.ast = 0.01;              % area of stapes
pa.arw = 0.01;              % area of round window
pa.ahe= 0;                  % area of helicotrema
% chamber sizes
pa.chsz=[1.312506 0.289117 0.001000 0.397377]; % err=28.14 28.14
pa.chpd=[1.000 -1.489 -0.589 -1.107 ; -0.408 1.000 0.251 -0.809 ; -0.611 -0.587 1.000 -0.859]; % err=28.08 28.02
% ------ parfit ------
pa.k1o=4.856e+08;
pa.r1o=5.648;
pa.m1o=0.0006127;
pa.k2o=1.718e+09;
pa.r2o=23.38;
pa.m2o=0.008588;
pa.k3o=3.085e+05;
pa.r3o=6.498;
pa.m3o=3.707;
pa.k4o=2.591e+05;
pa.r4o=0.02468;
pa.aco=0.007308;
pa.k1e=-4.2447;
pa.r1e=-0.7605;
pa.m1e=-0.0815;
pa.k2e=-4.3598;
pa.r2e=-0.1125;
pa.m2e=-0.0119;
pa.k3e=-0.0562;
pa.r3e=-0.0166;
pa.m3e=-0.0006;
pa.k4e=-4.3972;
pa.r4e=0.0001;
pa.ace=-0.4818;
pa.k1q=0.000004;
pa.r1q=0.000042;
pa.m1q=0.000006;
pa.k2q=0.000003;
pa.r2q=0.000086;
pa.m2q=0.000009;
pa.k3q=0.000018;
pa.r3q=-0.000040;
pa.m3q=0.000024;
pa.k4q=0.000080;
pa.r4q=0.000015;
pa.acq=0.000023;
pa.hbt=-12.0672;
% ------
end

%==============================================================

function [xx,Yb,Pd,Db,Dh,Hh]=modset(flst,pa)
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
% useful constructs
s=2i*pi*f;
dx=xl/(n-1);
srd=s*rho*dx;
% compute admittance for all x
x=transpose(linspace(0,xl,n));
ac=pa.aco*exp(pa.ace*x+pa.acq*x.^2); % area taper
bw=pa.bwo*exp(pa.bwe*x);             % BM width
[z1,z2,zo,za]=imped(x,s,pa);
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
