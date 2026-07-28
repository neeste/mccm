% parfit23 - fit 2 or 3 chamber model to MAP + forward latency
function parfit23
nc=3; fprintf('%d chambers\n',nc);
pa=modpar23(nc); % fetch model parameters
flst=500*2.^(-1:5); % frequency list
[xpk,mpk,ppk]=pk_tgt(flst);
%pa=fitpar23(pa,flst,xpk,mpk,ppk,2); % apex
%pa=fitpar23(pa,flst,xpk,mpk,ppk,1); % base
for k=1:20
    % fit model parameters to data
    pa=fitpar23(pa,flst,xpk,mpk,ppk,0);
    % compute frequency-domain cochlear model
    [xx,Yb,Pd,Db,Dh,Zc,Gf]=fdmod23(pa,flst);
    % plot model results
    fdm_plt(flst,xx,Yb,Pd,Db,Dh,pa.hbt,xpk,mpk,ppk);
    midear(8,flst,pa,Zc,Gf);
end
end

%------------------------------------------------

function fdm_plt(flst,xx,Yb,Pd,Db,Dh,hbt,xpk,mpk,ppk)
Zb=1./Yb;
Hh=Dh./Db;
mpxplt(1,xx,Yb,[ -60   40],[-0.3   0.3],'BM admittance');
rixplt(2,xx,Zb,[-0.5  1.5],[-1.5   0.5],'BM impedance');
mpxplt(3,xx,Pd,[ -30   70],[-26    1  ],'BM pressure difference');
mpxplt(4,xx,Db,[ -40   60],[-26    1  ],'BM displacement (nm at 0 dB SPL)');
mpxplt(5,xx,Dh,[ -40   60],[-26    1  ],'HB displacement (nm at 0 dB SPL)');
mpxplt(6,xx,Hh,[ -35   15],[ -0.1  0.5],'HB/BM filter (dB)');
excpat(7,xx,Dh,[-100    0],[-26    1  ],'excitation pattern',hbt,flst,xpk,mpk,ppk);
drawnow
end

function mpxplt(fig,xx,yy,dblim,phlim,lab)
mg=20*log10(abs(yy));
ph=angle(yy);
[nx,nf]=size(yy);
for k=1:nf
    phk=unwrap(ph(:,k))/(2*pi);
    phr=max(phk)-min(phk);
    if (phr<1.01)
        phc=round((max(phk)+min(phk))/2);
        phk=phk-phc;
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

function excpat(fig,xx,Dh,dblim,phlim,lab,hbt,~,xpk,mpk,ppk)
tpi=2*pi;
nx=length(xx);
db1=20*log10(abs(Dh))-hbt;
ph1=unwrap(angle(Dh))/tpi;
ph1(db1>120)=nan;
ph1=ph1-repmat(round(ph1(2,:)),nx,1);
% find peaks
hbxpk=xpk;
hbmpk=mpk;
hbppk=ppk;
for k=1:length(xpk)
    [~,ix]=max(db1(:,k));
    hbxpk(k)=(ix-1)/(nx-1);
    hbmpk(k)=db1(ix,k);
    hbppk(k)=ph1(ix,k);
end
% plot
figure(fig);clf;
subplot(211);
plot(xx,db1,xpk,mpk,'o',hbxpk,hbmpk,'k.')
ax=gca;ax.ColorOrderIndex=1;
axis([0 1 dblim])
ylabel('magnitude (dB)')
title(lab);
subplot(212)
plot(xx,ph1,xpk,ppk,'o',hbxpk,hbppk,'k.')
axis([0 1 phlim])
xlabel('x/L')
ylabel('phase (cyc)')
%write_data('moh24px.txt',[xpk' mpk' ppk']);
%write_data('moh24ep.txt',[xx db1 ph1]);
end

%------------------------------------------------

function [xpk,mpk,ppk]=pk_tgt(flst)
xpk= [ 0.80  0.72  0.59  0.47  0.33  0.20  0.05];
mpk=-[18.20  9.70  8.80 15.00 12.30 19.10 59.10];
ppk=-[10.55 15.52 17.96 20.02 22.08 21.83 13.80];
fpk=500*2.^(-1:5);
xpk=interp1(fpk,xpk,flst);
mpk=interp1(fpk,mpk,flst);
ppk=interp1(fpk,ppk,flst);
%xpk=xpk+0.05; % add hook ???
end

%------------------------------------------------

function pa=fitpar23(pa,flst,xpk,mpk,ppk,ot)
pasv=pa;
mxev=120;
if (ot==0)
    pv=getpar(pa);
    mxev=960;
elseif (ot==1)
    pv(1)=pa.kme;
    pv(2)=pa.rme;
elseif (ot==2)
    pv(1)=pa.khe;
    pv(2)=pa.rhe;
end
% fit parameters to neural data
op=optimset('MaxFunEvals',mxev);
ipk=1+round(xpk*(pa.n-1));
e1=refdev(pv,pa,flst,ipk,mpk,ppk,ot);
pv=fminsearch(@(pv) refdev(pv,pa,flst,ipk,mpk,ppk,ot),pv,op);
e2=refdev(pv,pa,flst,ipk,mpk,ppk,ot);
pa=pasv;
if (ot==0)
    pa=setpar(pa,pv);
    prnpar(pa);
elseif (ot==1)
    pa.kme=pv(1);
    pa.rme=pv(2);
    fprintf('pa.kme=%.4g; pa.rme=%.2f; pa.mme=%.4g; %% ',pv,pa.mme)
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
    no=10; ne=10; nq=10;
    io=1:no; ie=(1:ne)+no; iq=(1:nq)+(no+ne);  % pv indices
    mnpv=min(pv(io));
    if (mnpv<2e-6) err=1e9; return; end % minimum mass
    err=err+mean(abs(pv(ie)))*8;        % minimize linear exponents
    err=err+mean(abs(pv(iq)))*8;        % minimize quadratic exponents
    pa=setpar(pa,pv);
    if (pa.gpo>1) err=err+pa.gpo*100; end % limit partition gain < 1
    mxex=12;
    lr1=log(pa.r1o) + pa.r1e + pa.r1q;
    lm1=log(pa.m1o) + pa.m1e + pa.m1q;
    lrm1=lr1-lm1;
    if (lrm1>mxex) err=err+(lrm1-mxex); end % limit tc1
    lr2=log(pa.r2o) + pa.r2e + pa.r2q;
    lm2=log(pa.m2o) + pa.m2e + pa.m2q;
    lrm2=lr2-lm2;
    if (lrm2>mxex) err=err+(lrm2-mxex); end % limit tc2
elseif (ot==1)
    if (min(pv(1:2))<0.0005) err=1e9; return; end
    if (pv(1)>1e6)           err=1e9; return; end
    if (pv(2)>1800)          err=1e9; return; end
    pa.kme=pv(1);
    pa.rme=pv(2);
elseif (ot==2)
    if (min(pv(1:2))<1e-6) err=1e9; return; end
    if (max(pv(1:2))>1e4)  err=1e9; return; end
    pa.khe=pv(1);
    pa.rhe=pv(2);
end
% compute fdm
[~,~,~,~,Dh]=fdmod23(pa,flst);
D1=20*log10(abs(Dh))-pa.hbt;
[nx,nf]=size(D1);
dx=round(nx*0.02);
shlst=[0.3 0.3 4 6 6 3 2];
sh=interp1(500*2.^(-1:5),shlst,flst);
err=err*nf;
for k=1:nf
    ii=D1(:,k)>-950;
    if (sum(ii)==0) err=1e6; return; end
    [~,ix]=max(D1(:,k));
    mx=D1(ix,k);
    ph=unwrap(angle(Dh(:,k)))/(2*pi);
    err=err+sum(ph>ph(1));                  % positive phase slope ?
    err=err+abs(mx-mpk(k))*2;               % peak
    err=err+abs(ph(ix)-ppk(k))*2;           % phase
    err=err+abs((ix-ipk(k)))*20;            % fp-map
    if (1<(ix-dx))
        ls=mx-D1(ix-dx,k);
        err=err+abs(ls-sh(k))*2;            % left shoulder
    end
    if ((ix+dx)<nx)
        rs=mx-D1(ix+dx,k);
        err=err+abs(rs-sh(k))*2;            % right shoulder
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
k=k+1;pv(k)=pa.k4o;
k=k+1;pv(k)=pa.aco;
k=k+1;pv(k)=pa.k1e;
k=k+1;pv(k)=pa.r1e;
k=k+1;pv(k)=pa.m1e;
k=k+1;pv(k)=pa.k2e;
k=k+1;pv(k)=pa.r2e;
k=k+1;pv(k)=pa.m2e;
k=k+1;pv(k)=pa.k3e;
k=k+1;pv(k)=pa.r3e;
k=k+1;pv(k)=pa.k4e;
k=k+1;pv(k)=pa.ace;
k=k+1;pv(k)=pa.k1q;
k=k+1;pv(k)=pa.r1q;
k=k+1;pv(k)=pa.m1q;
k=k+1;pv(k)=pa.k2q;
k=k+1;pv(k)=pa.r2q;
k=k+1;pv(k)=pa.m2q;
k=k+1;pv(k)=pa.k3q;
k=k+1;pv(k)=pa.r3q;
k=k+1;pv(k)=pa.k4q;
k=k+1;pv(k)=pa.acq;
k=k+1;pv(k)=pa.gpo;
end

function pa=setpar(pa,pv)
k=0;
k=k+1;pa.k1o=exp(round(log(pv(k)),6));
k=k+1;pa.r1o=exp(round(log(pv(k)),6));
k=k+1;pa.m1o=exp(round(log(pv(k)),6));
k=k+1;pa.k2o=exp(round(log(pv(k)),6));
k=k+1;pa.r2o=exp(round(log(pv(k)),6));
k=k+1;pa.m2o=exp(round(log(pv(k)),6));
k=k+1;pa.k3o=exp(round(log(pv(k)),6));
k=k+1;pa.r3o=exp(round(log(pv(k)),6));
k=k+1;pa.k4o=exp(round(log(pv(k)),6));
k=k+1;pa.aco=exp(round(log(pv(k)),6));
k=k+1;pa.k1e=round(pv(k),4);
k=k+1;pa.r1e=round(pv(k),4);
k=k+1;pa.m1e=round(pv(k),4);
k=k+1;pa.k2e=round(pv(k),4);
k=k+1;pa.r2e=round(pv(k),4);
k=k+1;pa.m2e=round(pv(k),4);
k=k+1;pa.k3e=round(pv(k),4);
k=k+1;pa.r3e=round(pv(k),4);
k=k+1;pa.k4e=round(pv(k),4);
k=k+1;pa.ace=round(pv(k),4);
k=k+1;pa.k1q=round(pv(k),6);
k=k+1;pa.r1q=round(pv(k),6);
k=k+1;pa.m1q=round(pv(k),6);
k=k+1;pa.k2q=round(pv(k),6);
k=k+1;pa.r2q=round(pv(k),6);
k=k+1;pa.m2q=round(pv(k),6);
k=k+1;pa.k3q=round(pv(k),6);
k=k+1;pa.r3q=round(pv(k),6);
k=k+1;pa.k4q=round(pv(k),6);
k=k+1;pa.acq=round(pv(k),6);
k=k+1;pa.gpo=round(pv(k),4);
end

function prnpar(pa)
fp=fopen('fdm23par.txt','wt');
fprintf(fp,'pa.k1o=%.6g;\n',pa.k1o);
fprintf(fp,'pa.r1o=%.6g;\n',pa.r1o);
fprintf(fp,'pa.m1o=%.6g;\n',pa.m1o);
fprintf(fp,'pa.k2o=%.6g;\n',pa.k2o);
fprintf(fp,'pa.r2o=%.6g;\n',pa.r2o);
fprintf(fp,'pa.m2o=%.6g;\n',pa.m2o);
fprintf(fp,'pa.k3o=%.6g;\n',pa.k3o);
fprintf(fp,'pa.r3o=%.6g;\n',pa.r3o);
fprintf(fp,'pa.k4o=%.6g;\n',pa.k4o);
fprintf(fp,'pa.r4o=%.6g;\n',pa.r4o);
fprintf(fp,'pa.aco=%.6g;\n',pa.aco);
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
fprintf(fp,'pa.k1q=%.6f;\n',pa.k1q);
fprintf(fp,'pa.r1q=%.6f;\n',pa.r1q);
fprintf(fp,'pa.m1q=%.6f;\n',pa.m1q);
fprintf(fp,'pa.k2q=%.6f;\n',pa.k2q);
fprintf(fp,'pa.r2q=%.6f;\n',pa.r2q);
fprintf(fp,'pa.m2q=%.6f;\n',pa.m2q);
fprintf(fp,'pa.k3q=%.6f;\n',pa.k3q);
fprintf(fp,'pa.r3q=%.6f;\n',pa.r3q);
fprintf(fp,'pa.k4q=%.6f;\n',pa.k4q);
fprintf(fp,'pa.r4q=%.6f;\n',pa.r4q);
fprintf(fp,'pa.acq=%.6f;\n',pa.acq);
fprintf(fp,'pa.hbt=%.4f;\n',pa.hbt);
fprintf(fp,'pa.gpo=%.4f;\n',pa.gpo);
fclose(fp);
end

%==============================================================

function [xx,Yb,Pd,Db,Dh,Zc,Gf]=fdmod23(pa,flst)
L=pa.xl;
N=pa.n;
nf=length(flst);
dx=L/(N-1);
nc=length(pa.chsz);
Yb=zeros(N,nf);
Pd=zeros(N,nf);
Db=zeros(N,nf);
Dh=zeros(N,nf);
Zc=zeros(1,nf);
% loop over frequencies
for k=1:nf
    f=flst(k);
    s=2i*pi*f;
    [a,q,x,y,D]=mxfill(f,pa);    % initialize arrays
    p=mxsolve(a,q);
    [pd,yb,dd]=p_dif(p,y,f,D);
    [pd,yb,dd]=p_nan(pd,yb,dd);  % restrict plotting range
    Pd(:,k)=pd;
    Yb(:,k)=yb;
    Db(:,k)=dd(:,1);
    Dh(:,k)=dd(:,2);
    Vs=((p(nc,2)-p(nc,1))/dx)/s; % stapes velocity
    Zc(k)=p(nc,1)/Vs;            % cochlear impedance
end
xx=x/pa.xl;
Zc=Zc/pa.ast; % convert from mechanical to acoustic impedance ???
Gf=Pd(1,:);
end

%--------------------------------------------------------------

function [a,q,x,y,D]=mxfill(f,pa)
% local parameter values
chsz=pa.chsz;
chsz=chsz*2/sum(chsz); % normalize channel sizes;
n=pa.n;
xl=pa.xl; rho=pa.rho; ast=pa.ast; g=pa.gam;
bwo=pa.bwo; bwe=pa.bwe;
gh=pa.gpo*ones(n,1);
% useful constructs
s=2i*pi*f;
dx=xl/(n-1);
srd=s*rho*dx;
% compute admittance for all x
x=transpose(linspace(0,xl,n));
[z1,z2,zh,za,ac]=imped(x,s,pa);
bw=bwo*exp(bwe*x);
m=length(pa.chsz);
A=zeros(n,m);
Y=zeros(n,m,m);
a=zeros(3,m,m,n);
y=zeros(n,m,m);
% initialize base & apex
Zm=pa.kme/s+pa.rme+pa.mme*s;
Zh=pa.khe/s+pa.rhe+pa.mhe*s;
z1(n)=Zh;z2(n)=1e-6;zh(n)=1e-6;za(n)=0;
% initialize partition admittance
zg=g*za;
zk=zeros(m,m);
if (m==2)
    D=[1 -1;0 0];                   % select PD
    B=[-1;1];                       % basal BC
elseif (m==3)
    D=[1 0 -1;0 [1 -1]; 0 0 0]; % select PD
    B=[-1;0;1];                 % basal BC
end
for k=1:n
    A(k,:)=ac(k)*chsz;
    if (m==2)
        hh=gh(k)*z2(k)/(z2(k)+zh(k));
        zk(1,1)=z1(k)+(gh(k)*zh(k)-zg(k))*hh;
    elseif (m==3)
        zk(1,1:2)=[z1(k)  gh(k)*zh(k)-zg(k)]; 
        zk(2,1:2)=[z2(k)  -(zh(k)+z2(k))];
    end
    zk(m,:)=1; % conserve fluid volume
    Y(k,:,:)=zk\D;
    if (m==2)
        y(k,1,1)=1/zk(1,1);   % select Ybm 
        y(k,2,:)=hh*y(k,1,:); % select Yhb
    elseif (m==3)
        y(k,:,:)=inv(zk);     % select Ybm & Yhb
    end
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
% fill matrix
A1=2*diag(A(1,:));
A2=-A1;
An=diag(A(n,:));
a(1,:,:,1)=zeros(m,m);
a(2,:,:,1)=A1;
a(3,:,:,1)=A2;
a(1,:,:,n)=-An;
a(2,:,:,n)=An;
a(3,:,:,n)=zeros(m,m);
% middle-ear boundary constraints
if (0)
    [~,~,A3]=midear(0,f,pa);
else
    A3=-2*(srd*ast./Zm);
end
a(2,:,:,1)=A1;
a(3,:,:,1)=A2;
% initialize q
q=zeros(m,n);
q(:,1)=A3*B;
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

function [pd,yd,dd]=p_dif(p,y,f,D)
pref=2.848e-4; % dB SPL pressure reference (dyn/cm^2)
dref=1e-7;     % 1-nm displacement reference (cm)
n=size(p,2);
pd=zeros(n,1);
yd=zeros(n,1);
vd=zeros(n,2);
for k=1:n
   pk=D*p(:,k);
   yk=squeeze(y(k,:,:));
   vk=yk*pk;
   vd(k,:)=vk(1:2)*(pref/dref); % BM,HB velocity
   pd(k)=pk(1);                 % ST-SV pressure difference
   yd(k)=vk(1)/pd(k);           % BM admittance
end
dd=vd./(2i*pi*f);               % BM,HB displacement
end

function [pd,yd,dd]=p_nan(pd,yd,dd)
% use nan to restrict plotted range
n=length(pd);
for k=1:n
   if(abs(pd(k,1))<0.005)   
      ii=k:n;
      pd(ii,:)=nan;
      yd(ii,:)=nan;
      dd(ii,:)=nan;
      break;
   end
end
end

%==============================================================

function pa=modpar23(nc)
% chamber sizes
if (nc==2)
    pa=modpar22c2;
elseif (nc==3)
    pa=modpar23c3;
end
pa.man=(pa.xl>3);
end

%--------------------------------------------------------------

function pa=modpar22c2
pa.chsz=[1 1];                 % two-chamber sizes
%----------------------------
pa.gam = 1;                    % NDR multiplier
pa.m = 1;                      % number of points across fluid
pa.n = 1401;                   % number of points along BM
pa.xl = 3.5;                   % scala length
pa.yw = 0.1;                   % scala width
pa.zh = 0.1;                   % scala height
pa.rho = 1;                    % fluid density
pa.bwo = 0.05;                 % BM width at base
pa.bwe = 0;                    % BM width taper
pa.gpo = 1;                    % partition gain (HB re BM)
pa.isv = [562 486 408 325 235 138 34]; % BM locations to save
pa.xtap=0.0; pa.xtex=6;
pa.hbt=0; pa.xp=1;
pa.khe=0.0001; pa.rhe=0.0001; pa.mhe=0.01;
pa.kme=0.01058; pa.rme=728.00; pa.mme=0.04; % err=796.4 795.9
% middle-ear parameters
pa.mco=30; pa.rco=1.2e6; pa.rrw=2e5; pa.krw=5e7;
pa.mma=0.017; pa.rma=80; pa.kma=3e5; pa.aed=0.33;
pa.rim=400; pa.kim=5e6; pa.gm=1;
pa.mst=0.017; pa.rst=80; pa.kst=3e5; pa.ast=0.01;
pa.mrw=5e-3; pa.rrw=20; pa.krw=1.5e5; pa.arw=0.0625;
% ---- prft ----
pa.k1o=1.05e+08;
pa.r1o=126.6;
pa.m1o=0.002144;
pa.k2o=9.91e+06;
pa.r2o=14.94;
pa.m2o=0.002031;
pa.k3o=1.878e+07;
pa.r3o=1.024;
pa.k4o=2.352e+07;
pa.r4o=0;
pa.aco=0.01;
pa.k1e=-1.0118;
pa.r1e=0.9895;
pa.m1e=1.7590;
pa.k2e=-2.2567;
pa.r2e=1.1982;
pa.m2e=-0.2543;
pa.k3e=-3.6542;
pa.r3e=2.5617;
pa.k4e=0.1840;
pa.r4e=0.0000;
pa.ace=0.0000;
pa.k1q=-0.281818;
pa.r1q=-0.359353;
pa.m1q=-0.493720;
pa.k2q=-0.457496;
pa.r2q=-0.615777;
pa.m2q=0.637336;
pa.k3q=0.770590;
pa.r3q=-0.385747;
pa.k4q=-0.568728;
pa.r4q=0.000000;
pa.acq=0.000000;
end

%--------------------------------------------------------------

function pa=modpar23c3
pa.chsz=[0.95 0.05 1];         % three-chamber sizes
%----------------------------
pa.gam = 1;                    % NDR multiplier
pa.m = 1;                      % number of points across fluid
pa.n = 701;                    % number of points along BM
pa.xl = 3.5;                   % scala length
pa.yw = 0.1;                   % scala width
pa.zh = 0.1;                   % scala height
pa.rho = 1;                    % fluid density
pa.bwo = 0.05;                 % BM width at base
pa.bwe = 0;                    % BM width taper
pa.gpo = 1;                    % partition gain (HB re BM)
pa.isv = [562 486 408 325 235 138 34]; % BM locations to save
pa.xtap=0; pa.xtex=6;
pa.hbt=60; pa.xp=1;
pa.khe=0.0001; pa.rhe=0.0001; pa.mhe=0.01;
pa.kme=4.878e+04; pa.rme=29.51; pa.mme=0.001; % err=91.44 91.25
% middle-ear parameters
pa.mco=30; pa.rco=1.2e6; pa.rrw=2e5; pa.krw=5e7;
pa.mma=0.017; pa.rma=80; pa.kma=3e5; pa.aed=0.33;
pa.rim=400; pa.kim=5e6; pa.gm=1;
pa.mst=0.017; pa.rst=80; pa.kst=3e5; pa.ast=0.01;
pa.mrw=0.005; pa.rrw=20; pa.krw=1.5e5; pa.arw=0.0625;
% ---- parfit ----
pa.k1o=454294;
pa.r1o=0.0357451;
pa.m1o=4.73134e-06;
pa.k2o=2210.34;
pa.r2o=91.9516;
pa.m2o=0.000374226;
pa.k3o=685320;
pa.r3o=0.000471628;
pa.k4o=743711;
pa.r4o=0;
pa.aco=0.106088;
pa.k1e=-2.6875;
pa.r1e=3.7457;
pa.m1e=-1.4550;
pa.k2e=-1.3110;
pa.r2e=-4.8000;
pa.m2e=-1.6839;
pa.k3e=-3.2455;
pa.r3e=0.6154;
pa.k4e=-5.2092;
pa.r4e=0.0000;
pa.ace=0.0133;
pa.k1q=0.183633;
pa.r1q=-1.600415;
pa.m1q=-0.617439;
pa.k2q=0.576781;
pa.r2q=1.652350;
pa.m2q=0.930405;
pa.k3q=-0.443910;
pa.r3q=0.205098;
pa.k4q=0.863907;
pa.r4q=0.000000;
pa.acq=0.003402;
pa.hbt=60.0000;
pa.gpo=0.9311;
end
