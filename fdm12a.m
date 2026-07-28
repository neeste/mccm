% fdm12a - Finite-difference solution of the 1-D cochlear model

function fdm12a
global n dx ast ase rho bw As g xl wd npk sbc rv icp
%
n=2501; 
nf=600;						% number of frequencies
g=1;							% NDR depth
rr=5e-5;						% roughness
utf=0.0;                % untapered fraction of scala area
wd=1;							% write data ?
%
% set parameter values
xl=3.5; yh=0.1; bw=0.1; 
rho=1; ast=0.01; ase=-5/xl; % ase=-1.43; % for man
x=transpose(linspace(0,xl,n));
dx=x(2)-x(1);
As=ast*(exp(ase*x)*(1-utf)+utf);
rand('state',0);
rv1=rr*rand(size(x));   % random variation in stiffness
rv0=zeros(size(rv1));   % no random variation in stiffness
sbc=0;                  % pressure=1 at x=0
%
% plot impedance, admittance, pressure, velocity & propagation
L=xl;rv=rv0;
xlim=[0 L];
zlim=[-3 1]*8000;
ylim=[-120 -60];
clim=[-7 1];
klim=[-60 120];
figure(1);clf
subplot(211);hold on;axis([xlim zlim])
plot(xlim,[0 0],'c:');
title('impedance');
ylabel('real');
subplot(212);hold on;axis([xlim zlim])
plot(xlim,[0 0],'c:');
xlabel('distance from stapes (cm)')
ylabel('imaginary');
figure(2);clf
subplot(211);hold on;axis([xlim ylim])
title('admittance');
ylabel('magnitude (dB)');
subplot(212);hold on;axis([xlim -0.5 0.5])
plot(xlim,[0.25 0.25],'c:');
xlabel('distance from stapes (cm)')
ylabel('phase (cyc)');
figure(3);clf
subplot(211);hold on;axis([xlim -40 60])
title('pressure ratio');
ylabel('level (dB)');
subplot(212);hold on;axis([xlim clim])
xlabel('distance from stapes (cm)')
ylabel('phase (cyc)')
figure(4);clf
subplot(211);hold on;axis([xlim -15 85])
title('velocity ratio')
ylabel('level (dB)')
subplot(212);hold on;axis([xlim clim])
xlabel('distance from stapes (cm)');
ylabel('phase (cyc)')
figure(5);clf
subplot(211);hold on;axis([xlim klim])
plot(xlim,[0 0],'c:')
title('propagation function')
ylabel('real')
subplot(212);hold on;axis([xlim klim])
xlabel('distance from stapes (cm)')
ylabel('imaginary')
drawnow
if (xl<3) o=1:6; else o=0:5; end
icp=zeros(1,length(o));
flst=250*(2.^o);
for j=1:length(flst)
   s=i*2*pi*flst(j);
   Zf=rho*s./As;			% scala impedance
   Yb=Ybm(x,s,g);			% BM admittance
   Zb=1./Yb;				% BM impedance
   P=press(Zf,Yb,x,0);	% scala pressure
   Vs=-((P(2)-P(1))/(x(2)-x(1))/Zf(1))/As(1);
   Vb=(Yb.*P/bw)/Vs;		% BM velocity
   kap=sqrt(Zf./Yb).*Yb;
   %
   Ydb=dbmag(Yb);
   Yph=phase(Yb);
   Pdb=dbmag(P/P(1));
   Pph=phase(P/P(1));
   Vdb=dbmag(Vb);
   Vph=phase(Vb);
   kr=real(kap);
   ki=imag(kap);
   figure(1)
   subplot(211);plot(x,real(Zb))
   subplot(212);plot(x,imag(Zb))
   figure(2)
   subplot(211);plot(x,Ydb)
   subplot(212);plot(x,Yph)
   figure(3)
   subplot(211);plot(x,Pdb)
   subplot(212);plot(x,Pph)
   figure(4)
   subplot(211);plot(x,Vdb)
   subplot(212);plot(x,Vph)
   figure(5)
   subplot(211);plot(x,kr)
   subplot(212);plot(x,ki)
   drawnow
   [xmx,imx]=max(Vdb);
   icp(j)=imx;
   write_data(sprintf('velo%d.txt',j),[x Vdb Vph Ydb Yph kr ki]);
end
fig=5;
%
% frequency-place map
fig=fig+1;figure(fig);clf;
xcp=x(icp);
f_cat=fpmap(x,0.456,0.021,0.80);
f_man=fpmap(x,0.165,0.021,0.99);
semilogy(xcp,flst,'o-',x,f_cat,x,f_man)
flim=[.05 20]*1000;
axis([xlim flim])
xlabel('distance from stapes (cm)');
ylabel('characteristic frequency (Hz)');
title('frequency-place map');
legend('model','cat','man')
%
% reflectance vs frequency
fig=fig+1;figure(fig);clf
flim=[50 20000];
o=8*(0:(nf-1))/nf-3;
f=1000*2.^o';
[rm1,rd1,pp1,pd1,vf1,kf]=refl(x,f,rv0,0,1);
[rm2,rd2,pp2,pd2,vf2,kf]=refl(x,f,rv1,0,2);
rm1=20*log10(rm1/100);
rm2=20*log10(rm2/100);
subplot(211);semilogx(f,rm1,f,rm2)
axis([flim -65 -15]);
title('reflectance');
ylabel('magnitude (dB)');
rd1(rd1<eps)=eps;
rd2(rd2<eps)=eps;
subplot(212);loglog(f,rd1,f,rd2)
axis([flim 1 40])
ylabel('delay (cyc)')
xlabel('frequency (Hz)')
drawnow
write_data('vf1.txt',[f dbmag(vf1) phase(vf1)]);
write_data('vf2.txt',[f dbmag(vf2) phase(vf2)]);
write_data('kf.txt', [f real(kf)   imag(kf)]);
%
% Qerb
fig=fig+1;figure(fig);clf
for k=1:length(icp)
   Vdb2=dbmag(vf2(:,k));
   [qe,bf]=Qerb(f,Vdb2);
   qq(k)=qe;
   ff(k)=bf;
end
f=1000*2.^(-6:0.1:6);
Qcat = 5.0*(f/1000).^0.37;	% from Shera et al. (2002)
Qman1=9.26*f./(f+230);		% from Shera et al. (2002)
Qman2=12.7*(f/1000).^0.30;	% from Shera et al. (2002)
loglog(ff,qq,'o-',f,Qcat,f,Qman1,f,Qman2)
axis([flim 1 30])
title('Qerb vs BF')
ylabel('quality factor')
xlabel('frequency (Hz)')
%legend('model','cat','man1','man2',4)
legend('model','cat','man1','man2')
drawnow
%
return

function [rm,rd,pp,pd,vf,kf]=refl(x,flst,rvx,ap,sn)
global rho bw As g rv dx icp
ap=1;
rv=rvx;  				% random variation in stiffness
r=zeros(size(flst));
pp=zeros(size(flst));
pph=zeros(size(flst));
vf=zeros(length(flst),length(icp));
kf=zeros(length(flst),length(icp));
for j=1:length(flst)
   s=i*2*pi*flst(j);
   Zf=rho*s./As;		% scala impedance
   Yb=Ybm(x,s,g);		% BM admittance
   [P,Pp,Pm]=press(Zf,Yb,x,ap);	% scala pressure
   r(j)=Pm(1)/Pp(1);
   pp(j)=max(dbmag(P/P(2)));		% peak pressure
   pph(j)=Pm(2);
   Vb=(Yb.*P/bw)/As(1);				% BM velocity
   kap=sqrt(Zf./Yb).*Yb;			% propagation
   vf(j,:)=Vb(icp).';
   kf(j,:)=kap(icp).';
end
rm=abs(r)*100; 			% reflectance magnitude (%)
rp=phase(r);				% reflectance phase (cyc)
rd=-dd(rp,log(flst));	% reflectance delay (cyc)
ph=phase(pph);				% retrograde pressure phase (cyc)
pd=-dd(ph,log(flst));	% retrograde pressure delay (cyc)
write_data(sprintf('refl%d.txt',sn),[flst rm rp rd pp pd]);
return

function [p,pp,pm]=press(z,y,x,ap)
global As
us=As(1);	% stapes velocity
if ap
   [p,pp,pm]=press_ap(z,y,x,us);
else
   [p,pp,pm]=press_ex(z,y,x,us);
end   
return

function [p,pp,pm]=press_ex(z,y,x,us)
global As npk sbc
n=length(x);
o=ones([n 1]);
q=zeros([n 1]);
dx=x(2)-x(1);
% differential log-impedance
dlnZ=-dd(log(abs(z)),x)*dx;
% setup fluid equations
aa=2+z.*y*dx*dx;
ap=-o-dlnZ/2;
am=-o+dlnZ/2;
a=spdiags([am aa ap],-1:1,n,n);
% at stapes: P'=-As*Zf
a(1,1)=1;
a(1,2)=0;
q(1)=1;
% at helicotrema P'=-sqrt(Yb*Zf)*P;
a(n,n-1)=-1;
a(n,n)=1+sqrt(y(n)*z(n))*dx;
q(n)=0;
% solve matrix equation
p=a\q;
% decompose P into P+ and P-
dp=dd(p,x);
kap=y.* sqrt(z./y);
pp=( p - dp ./ kap ) / 2;
pm=( p + dp ./ kap ) / 2;
return

function [p,pp,pm]=press_ap(z,y,x,us)
dx=x(2)-x(1);
y0=sqrt(y./z);
eps=dd(log(y0),x)/2;
kap=y.* sqrt(z./y);
ckd=cumsum(2*kap.*dx);
eee=-eps.*exp(-ckd).*dx;
ref=exp(ckd).*(sum(eee)-cumsum(eee));
pp0=us/(y0(1)-y0(1)*ref(1));
pp=pp0.*exp(-cumsum((eps+kap).*dx));
pm=pp.*ref;
p=pp+pm;
return

function y=dbmag(x)
y=20*real(log10(max(eps,x)));
return

function y=phase(x)
y=unwrap(angle(x))/(2*pi);
return

function pp=dd(p,x)	% gradient
n=length(p);
dp=p;
dx=x;
dp(1)=p(2)-p(1);
dx(1)=x(2)-x(1);
dp(2:(n-1))=p(3:n)-p(1:(n-2));
dx(2:(n-1))=x(3:n)-x(1:(n-2));
dp(n)=p(n)-p(n-1);
dx(n)=x(n)-x(n-1);
pp=dp./dx;
return

% BM admittance for 1D model
function y=Ybm(x,s,g)
global xl
% compensate for scala area decrease
if (xl<3)   % cat
   y=0.2./Zbm(x,s,g*1.667);
else        % man
   g=g*1.42;
   x=x+0.8;
   %x=x+0.75+5.5e-4*x.^6; % apical bend
   y=0.1./Zbm(x,s,g);
end
return

% BM impedance from 2D model
function z=Zbm(x,s,g)
global rv ase
Ke=ase;
Re=0;
Me=-ase;
K1=1e9*exp(Ke*x).*(1+rv);
R1=450*exp(Re*x);
M1=0.003*exp(Me*x);
K2=5e2*exp(Ke*x);
R2=1.0004*exp(Re*x);
M2=6e-9*exp(Me*x);
R3=-1*exp(Re*x);
%
Z1=K1/s + R1 + M1*s;
Z2=K2/s + R2 + M2*s;
Z3=R3;
z=Z1+g.*Z2.*Z3./(Z2+Z3);
return

function cf=fpmap(x,a,b,c)
cf=1000*a*(10.^(100*b*(1-x/max(x)))-c);
return

function [Q,BF]=Qerb(f,G)
[Gmx,m]=max(G);
BF=f(m(1));
ERB=trapz(f,10.^((G-Gmx)/10));
Q=BF/ERB;
return

function write_data(name,data)
global wd
if (wd)
	fid=fopen(name,'wt');
	fprintf(fid,'; %s\n',name);
	[m,n]=size(data);
	for i=1:m
   	for j=1:n
      	fprintf(fid,'%10.4g ',data(i,j));
   	end
      fprintf(fid,'\n');
   end
	fclose(fid);
end
return
