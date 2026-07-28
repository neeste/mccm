% parfit22.m - fit 2-chamber model parameters to AN threshold
function parfit22
% compute 2DOF cochlear model
pa=modpar22;
flst=500*2.^(-1:5); % frequency list
%Dt=10.^(-tc_x(flst,pa.n,pa.xp,1)/20);
%Dt=adj_thr(flst,Dt);
[xpk,mpk,ppk]=pk_tgt;
%pa=fitpar22(pa,flst,xpk,mpk,ppk,2); % apex
%pa=fitpar22(pa,flst,xpk,mpk,ppk,1); % base
for k=1:200
    % fit model parameters to data
    pa=fitpar22(pa,flst,xpk,mpk,ppk,0);
    % compute frequency-domain cochlear model
    [xx,Yb,Pd,Db,Dh,Zc,Gf]=fdmod22(pa,flst);
    % plot model results
    fdm_plt(flst,xx,Yb,Pd,Db,Dh,pa.hbt,xpk,mpk,ppk);
    midear(8,flst,pa,Zc,Gf);
end
end

%------------------------------------------------

function fdm_plt(flst,xx,Yp,Pd,Db,Dh,hbt,xpk,mpk,ppk)
Zp=1./Yp;
Hh=Dh./Db;
mpxplt(1,xx,Yp,[-100  -20],[ -0.5  0.5],'BM admittance');
rixplt(2,xx,Zp,[-950  950],[-950  950 ],'BM impedance');
mpxplt(3,xx,Pd,[ -60   40],[ -7    1  ],'BM pressure difference');
mpxplt(4,xx,Db,[-100    0],[ -7    1  ],'BM displacement (nm at 0 dB SPL)');
mpxplt(5,xx,Dh,[-100    0],[ -7    1  ],'HB displacement (nm at 0 dB SPL)');
mpxplt(6,xx,Hh,[ -30   20],[ -0.3  0.3],'HB/BM filter (dB)');
excpat(7,xx,Dh,[-100    0],[-25    1  ],'excitation pattern',hbt,flst,xpk,mpk,ppk);
drawnow
end

function mpxplt(fig,xx,yy,dblim,phlim,lab)
nx=length(xx);
mg=20*log10(abs(yy));
ph=angle(yy);
[nx,nf]=size(yy);
for k=1:nf
    phk=unwrap(ph(:,k))/(2*pi);
    phr=max(phk)-min(phk);
    if (phr<1.01)
        phk=phk-(max(phk)+min(phk))/2;
    end
    ph(:,k)=phk;
end
ph=ph-repmat(round(ph(2,:)),nx,1);
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

function excpat(fig,xx,Dh,dblim,phlim,lab,hbt,flst,xpk,mpk,ppk)
tpi=2*pi;
nx=length(xx);
nf=length(flst);
db1=20*log10(abs(Dh))-hbt;
ph1=unwrap(angle(Dh))/tpi;
ph1(db1>120)=nan;
ph1=ph1-repmat(round(ph1(2,:)),nx,1);
figure(fig);clf;
subplot(211);
plot(xx,db1,xpk,mpk,'o')
ax=gca;ax.ColorOrderIndex=1;
axis([0 1 dblim])
ylabel('magnitude (dB)')
title(lab);
subplot(212)
plot(xx,ph1,xpk,ppk,'o')
axis([0 1 phlim])
xlabel('x/L')
ylabel('phase (cyc)')
text(0.1,-6,'d1-d2')
end

%------------------------------------------------

function [xpk,mpk,ppk]=pk_tgt
xpk= [ 0.80  0.72  0.59  0.47  0.33  0.20  0.05];
mpk=-[18.20  9.70  8.80 15.00 12.30 19.10 59.00];
ppk=-[ 2.50  3.80  5.80  8.70 13.00 20.00 24.00];
end

%------------------------------------------------

function pa=fitpar22(pa,flst,xpk,mpk,ppk,ot)
pasv=pa;
mxev=120;
if (ot==0)
    pv=getpar(pa);
    mxev=96000;
elseif (ot==1)
    pv(1)=pa.kme;
    pv(2)=pa.rme;
elseif (ot==2)
    pv(1)=pa.khe;
    pv(2)=pa.rhe;
end
% fit parameters to neural data
ipk=1+round(xpk*(pa.n-1));
e1=refdev(pv,pa,flst,ipk,mpk,ppk,ot);
op=optimset('MaxFunEvals',mxev);
%pv=fminsearch(@(pv) refdev(pv,pa,flst,ipk,mpk,ppk,ot),pv,op);
e2=refdev(pv,pa,flst,ipk,mpk,ppk,ot);
pa=pasv;
if (ot==0)
    pa=setpar(pa,pv);
    prnpar(pa);
elseif (ot==1)
    pa.kme=pv(1);
    pa.rme=pv(2);
    fprintf('pa.kme=%.4g; pa.rme=%.2f; pa.mme=%.2f; %% ',pv,pa.mme)
elseif (ot==2)
    pa.khe=pv(1);
    pa.rhe=pv(2);
    fprintf('pa.khe=%.4g; pa.rhe=%.4g; pa.mhe=%.4g; %% ',pv,pa.mhe)
end
fprintf('err=%.4g %.4g\n',e1,e2)
end

function err=refdev(pv,pa,flst,ipk,mpk,ppk,ot)
err=0;
if (ot==0)
    no=9; ne=9; nq=9;
    io=1:no; ie=(1:ne)+no; iq=(1:nq)+(no+ne);  % pv indices
    if (min(pv(io))<2e-3) err=1e9; return; end % minimum mass
    err=err+mean(abs(pv(ie)));          % minimize linear exponents
    err=err+mean(abs(pv(iq)));          % minimize quadratic exponents
    %err=err+log(pv(9)/1e6)*40;          % minimize k4o
    %err=err+sqrt(pv(1)*pv(3)/pv(2))/2;  % minimize Qbm
    %err=err+sqrt(pv(4)*pv(6)/pv(5))/1;  % minimize Qsf
    pa=setpar(pa,pv);
elseif (ot==1)
    if (min(pv(1:2))<0.01) err=1e9; return; end
    if (pv(1)>1)           err=1e9; return; end
    if (pv(2)>1800)        err=1e9; return; end
    pa.kme=pv(1);
    pa.rme=pv(2);
elseif (ot==2)
    if (min(pv(1:2))<1e-6) err=1e9; return; end
    if (max(pv(1:2))>1e4)  err=1e9; return; end
    pa.khe=pv(1);
    pa.rhe=pv(2);
end
% compute fdm22
[~,~,~,~,Dh]=fdmod22(pa,flst);
D1=20*log10(abs(Dh))-pa.hbt;
nf=length(flst);
for k=1:nf
    ii=D1(:,k)>-950;
    if (sum(ii)==0) err=1e6; return; end
    [mx,ix]=max(D1(:,k));
    ph=unwrap(angle(Dh(:,k)))/(2*pi);
    err=err+abs(mx-mpk(k));
    err=err+abs(ph(ix)-ppk(k))*80;
    err=err+abs((ix-ipk(k)))*10;
    %qqq=abs(min(diff(diff(D1(ii,k)))))*10; % sharpness of D1 peak
    %err=err+qqq;
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
k=k+1;pv(k)=pa.k4o;
k=k+1;pv(k)=pa.k1e;
k=k+1;pv(k)=pa.r1e;
k=k+1;pv(k)=pa.m1e;
k=k+1;pv(k)=pa.k2e;
k=k+1;pv(k)=pa.r2e;
k=k+1;pv(k)=pa.m2e;
k=k+1;pv(k)=pa.k3e;
k=k+1;pv(k)=pa.r3e;
k=k+1;pv(k)=pa.k4e;
k=k+1;pv(k)=pa.k1q;
k=k+1;pv(k)=pa.r1q;
k=k+1;pv(k)=pa.m1q;
k=k+1;pv(k)=pa.k2q;
k=k+1;pv(k)=pa.r2q;
k=k+1;pv(k)=pa.m2q;
k=k+1;pv(k)=pa.k3q;
k=k+1;pv(k)=pa.r3q;
k=k+1;pv(k)=pa.k4q;
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
k=k+1;pa.k4o=pv(k);
k=k+1;pa.k1e=pv(k);
k=k+1;pa.r1e=pv(k);
k=k+1;pa.m1e=pv(k);
k=k+1;pa.k2e=pv(k);
k=k+1;pa.r2e=pv(k);
k=k+1;pa.m2e=pv(k);
k=k+1;pa.k3e=pv(k);
k=k+1;pa.r3e=pv(k);
k=k+1;pa.k4e=pv(k);
k=k+1;pa.k1q=pv(k);
k=k+1;pa.r1q=pv(k);
k=k+1;pa.m1q=pv(k);
k=k+1;pa.k2q=pv(k);
k=k+1;pa.r2q=pv(k);
k=k+1;pa.m2q=pv(k);
k=k+1;pa.k3q=pv(k);
k=k+1;pa.r3q=pv(k);
k=k+1;pa.k4q=pv(k);
end

function prnpar(pa)
fp=fopen('fdm22par.txt','wt');
fprintf(fp,'pa.k1o=%.4g;\n',pa.k1o);
fprintf(fp,'pa.r1o=%.4g;\n',pa.r1o);
fprintf(fp,'pa.m1o=%.4g;\n',pa.m1o);
fprintf(fp,'pa.k2o=%.4g;\n',pa.k2o);
fprintf(fp,'pa.r2o=%.4g;\n',pa.r2o);
fprintf(fp,'pa.m2o=%.4g;\n',pa.m2o);
fprintf(fp,'pa.k3o=%.4g;\n',pa.k3o);
fprintf(fp,'pa.r3o=%.4g;\n',pa.r3o);
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
fprintf(fp,'pa.k4q=%.6f;\n',round(pa.k4q,6));
fprintf(fp,'pa.r4q=%.6f;\n',round(pa.r4q,6));
fprintf(fp,'pa.acq=%.6f;\n',round(pa.acq,6));
fclose(fp);
end

%================================================

function pa=modpar22
pa.gam = 1;                    % NDR multiplier
pa.m = 1;                      % number of points across fluid
pa.n = 1401;                   % number of points along BM
pa.xl = 3.5;                   % scala length
pa.yw = 0.1;                   % scala width
pa.zh = 0.1;                   % scala height
pa.rho = 1;                    % fluid density
pa.bwo = 0.05;                 % BM width at base
pa.bwe = 0;                    % BM width taper
pa.isv = [562 486 408 325 235 138 34]; % BM locations to save
pa.xtap=0.0; pa.xtex=6;
pa.hbt=0; pa.xp=1;
pa.khe=0.0001; pa.rhe=0.0001; pa.mhe=0.01; % err=23.21 23.21
pa.kme=0.01; pa.rme=1000; pa.mme=0.01; % err=23.21 23.21
% middle-ear parameters
pa.mco=30; pa.rco=1.2e6; pa.rrw=2e5; pa.krw=5e7;
pa.mma=0.017; pa.rma=80; pa.kma=3e5; pa.aed=0.33;
pa.rim=400; pa.kim=5e6; pa.gm=1;
pa.mst=0.017; pa.rst=80; pa.kst=3e5; pa.ast=0.01;
pa.mrw=5e-3; pa.rrw=20; pa.krw=1.5e5; pa.arw=0.0625;
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
Gf=Pd(1,:)/2;
return
end

function [a,q,x,Yb,Hh,D]=mxfill(f,pa)
% set parameter values
W=pa.yw;
H=pa.zh;
N=pa.n;
gm=pa.gam;
% useful constructs
s=2i*pi*f;
dx=pa.xl/(pa.n-1);
srd=s*pa.rho*dx;
% partition admittance
x=transpose(linspace(0,pa.xl,N));
[z1,z2,z3,z4,ac]=imped(x,s,pa);
Hh=z2./(z2+z3);
%Zm=pa.kme/s+pa.rme+pa.mme*s;
Zh=pa.khe/s+pa.rhe+pa.mhe*s;
Zb=z1+(z3-gm*z4).*Hh;
Zb(N-1)=Zh;
Yb=1./Zb;
% middle-ear boundary constraints
if (1)
    [A1,A2,A3]=midear(0,f,pa);
    A2=-A2/srd;
    A1=A1-A2;
    A3=A3*0.05;
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
%yr=(dx/H).*Yb;
yr=Yb.*(dx*W./ac);
a=[-o 2+srd*yr -o];
a(  1,1)=A1; a(1,2)=A2; % basal BC
a(N-1,N)=-1; a(N,N)=1;  % apical BC
a=spdiags(a,-1:1,N,N);
% initialize q
q=zeros(N,1);
q(1)=A3;
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
