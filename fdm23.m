% fdm23 - multi-chamber frequency-domain model of cochlea
function fdm23
nch=3;                    % number of channels
pa=modpar23(nch);         % fetch model parameters
flst1=500*2.^(-1:5);      % frequency list 1
flst2=500*2.^(-2:0.1:5);  % frequency list 2
[xpk,mpk,ppk]=pk_tgt(flst1);
report_err(pa,flst1,xpk,mpk,ppk);
[xx,Yb,Pd,Db,Dh]=fdmod23(pa,flst1);
% plot excitation patterns
fdm_plt(flst1,xx,Yb,Pd,Db,Dh,pa.hbt,xpk,mpk,ppk,nch);
% plot tuning curves
[xx,~,~,Db,Dh,Zc,Gf]=fdmod23(pa,flst2);
Ze=midear_imp(pa,flst2,Zc);
midear_plt(8,flst2,Zc,Gf,Ze);
dtc_plt( 9,xx,flst1,flst2,Db,pa,xpk,mpk,ppk,nch,'bm');
dtc_plt(10,xx,flst1,flst2,Dh,pa,xpk,mpk,nch,ppk,'hb');
end

%------------------------------------------------

function fdm_plt(flst,xx,Yp,Pd,Db,Dh,hbt,xpk,mpk,ppk,m)
if (m<1) return; end
Zp=1./Yp;
Hh=Dh./Db;
if (m<3)
    mpxplt(1,xx,Yp,[-100    0],[ -0.3  0.5],'BM admittance');
    rixplt(2,xx,Zp,[-200  200],[-200   200],'BM impedance');
    mpxplt(3,xx,Pd,[ -50   50],[-9     1  ],'BM pressure difference');
    mpxplt(4,xx,Db,[-100    0],[-9     1  ],'BM displacement (nm at 0 dB SPL)');
    mpxplt(5,xx,Dh,[-100    0],[-9     1  ],'HB displacement (nm at 0 dB SPL)');
    mpxplt(6,xx,Hh,[ -35   15],[-0.3   0.3],'HB/BM filter (dB)');
    excpat(7,xx,Dh,[-100    0],[-9     1  ],'excitation pattern',hbt,flst,xpk,mpk,ppk);
else
    mpxplt(1,xx,Yp,[ -60   40],[ -0.3  0.3],'BM admittance');
    rixplt(2,xx,Zp,[-0.5  1.5],[-1.5   0.5],'BM impedance');
    mpxplt(3,xx,Pd,[ -30   70],[-9     1  ],'BM pressure difference');
    mpxplt(4,xx,Db,[ -40   60],[-9     1  ],'BM displacement (nm at 0 dB SPL)');
    mpxplt(5,xx,Dh,[ -40   60],[-9     1  ],'HB displacement (nm at 0 dB SPL)');
    mpxplt(6,xx,Hh,[ -35   15],[-0.5   1  ],'HB/BM filter (dB)');
    excpat(7,xx,Dh,[-100    0],[-9     1  ],'excitation pattern',hbt,flst,xpk,mpk,ppk);
end
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
[nx,nf]=size(Dh);
db1=20*log10(abs(Dh))-hbt;
ph1=unwrap(angle(Dh))/tpi;
ph1(db1>120)=nan;
ph1=ph1-repmat(round(ph1(2,:)),nx,1);
% find peaks
kk=1:(nx-1);
hbxpk=zeros(1,nf);
hbmpk=zeros(1,nf);
hbppk=zeros(1,nf);
for k=1:length(xpk)
    [~,ix]=max(db1(kk,k));
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
ppk(7)=nan;
plot(xx,ph1,xpk,ppk,'o',hbxpk,hbppk,'k.')
axis([0 1 phlim])
xlabel('x/L')
ylabel('phase (cyc)')
write_data('moh24px.txt',[xpk' mpk' ppk']);
write_data('moh24ep.txt',[xx db1 ph1]);
end

%------------------------------------------------

function dtc_plt(fig,xx,flst1,flst2,Dh,pa,xpk,mpk,ppk,m,loc)
if (m<1) return; end
lab='displacement tuning curve';
tuncrv(fig,xx,flst1,flst2,Dh,[-100 20],[-9 1],lab,loc,pa,xpk,mpk,ppk);
end

function tuncrv(fig,xx,flst1,flst2,Dh,dblim,phlim,lab,loc,pa,xpk,mpk,ppk)
hbt=pa.hbt;
lab1=sprintf('nch=%d',pa.m);
lab2=pa.parlab;
flim=[0.2 20];
tpi=2*pi;
fk1=flst1/1000;
fk2=flst2/1000;
nx=length(xx);
ii=round(1+xpk*(nx-1));
nn=length(ii);
db1=20*log10(abs(Dh))-hbt;
ph1=unwrap(angle(Dh))/tpi;
ph1(db1>120)=nan;
ph1=ph1-repmat(round(ph1(2,:)),nx,1);
db2=db1(ii,:);
ph2=ph1(ii,:);
dbpk=ones(1,nn);
phpk=ones(1,nn);
fkpk=ones(1,nn);
for k=1:nn
    [~,kk]=max(db2(k,:));
    dbpk(k)=db2(k,kk);
    phpk(k)=ph2(k,kk);
    fkpk(k)=fk2(kk);
end
figure(fig);clf;
subplot(2,1,1);
semilogx(fk2,db2,fk1,mpk,'ro',fkpk,dbpk,'b.')
ax=gca;ax.ColorOrderIndex=1;
axis([flim dblim])
ylim(dblim)
ylabel('magnitude (dB)')
title([loc ' ' lab]);
text(10,0,lab2)
subplot(2,1,2)
ppk(7)=nan;
semilogx(fk2,ph2,fk1,ppk,'ro',fkpk,phpk,'b.')
axis([flim phlim])
ylim(phlim)
xlabel('frequency')
ylabel('phase (cyc)')
text(0.3,-20,lab1)
drawnow
fn1=sprintf('moh24%spf.txt',loc);
fn2=sprintf('moh24%stc.txt',loc);
write_data(fn1,[fk1' mpk' ppk' fkpk' dbpk' phpk']);
write_data(fn2,[fk2' db2' ph2']);
end

function midear_plt(fig,f,Zc,Gf,Ze)
if (fig==0) return; end
% Ze - middle-ear impedance (Hajicek 2024, ARO)
% Gf - middle-ear forward gain (Puria 2003, JASA)
% Zc - cochlea input impedance (Puria 2003, JASA)
ef=[0.1 0.2 0.5 1 2 5 10];                      % (kHz)
emZe=10^7*[22.8 13.3 6.5 4.8 4.1 4.0 4.9];      % (mks ohm)
emZc=1e9*[8.36 6.28 9.936 10.7 15.7 22.6 19.8]; % (mks ohm)
emGf=[-5.63 -2.77 11.0 16.6 11.7 0.25 -10.7];   % (dB)
%
fk=f/1000;
Zc=Zc*10^5;       % convert ftom cgs to mks
Ze=Ze*10^5;       % convert ftom cgs to mks
if (nargout<1)
    % plot
    figure(fig);clf
    subplot(2,2,1)
    loglog(fk,abs(Ze),ef,emZe);
    xlim([0.1 20])
    ylabel('Ze (mks ohm)')
    title('middle-ear impedance')
    subplot(2,2,2)
    loglog(fk,abs(Zc),ef,emZc);
    xlim([0.1 20])
    ylim([10^9 10^11])
    ylabel('Zc (mks ohm)')
    title('cochlea impedance')
    subplot(2,2,3)
    dbGf=20*log10(abs(Gf));
    semilogx(fk,dbGf,ef,emGf);
    xlim([0.1 20])
    ylim([-40 40])
    xlabel('frequency (Hz)')
    ylabel('Ps / Pe (dB)')
    title('middle-ear gain')
    drawnow
end
end

function Ze=midear_imp(pa,f,Zc)
s=2i*pi*f;
ast=pa.ast;
ama=pa.ama;
aed=pa.aed;
gme=pa.gme;
Zma=s*pa.mma+pa.rma+pa.kma./s;
Zim=pa.rim+pa.kim./s;
Zst=s*pa.mst+pa.rst+pa.kst./s;
Zed=s*pa.med+pa.red+pa.ked./s;
Zii=1./(1./Zim+1./(Zst+ast*Zc));
Zmi=Zma+gme*Zii;
Ze=1./(1./(ama*Zmi)+1./(aed*Zed));
end

%------------------------------------------------

function [xpk,mpk,ppk]=pk_tgt(flst)
xpk= [ 0.80  0.72  0.59  0.47  0.33  0.20  0.05];
mpk=-[18.20  9.70  8.80 15.00 12.30 19.10 59.10];
%ppk=-[10.55 15.52 17.96 20.02 22.08 21.83 13.80];
fpk=500*2.^(-1:5);
ppk=-min(4,12*(fpk/1000).^0.5);
xpk=interp1(fpk,xpk,flst);
mpk=interp1(fpk,mpk,flst);
ppk=interp1(fpk,ppk,flst);
%xpk=xpk+0.05; % add hook ?
end

%------------------------------------------------

function write_data(fn,data)
[nr,nc] = size(data);
fp=fopen(fn,'wt');
fprintf(fp,'; %s\n', fn);
for i=1:nr
    for j=1:nc
        fprintf(fp,' %14.5g',data(i,j));
    end
    fprintf(fp,'\n');
end
fclose(fp);
end

%------------------------------------------------

function report_err(pa,flst,xpk,mpk,ppk)
fprintf('nch=%d nme=%d\n',pa.m,pa.nmev);
if (pa.m<1) return; end
pv=getpar(pa);
ipk=1+round(xpk*(pa.n-1));
err=refdev(pv,pa,flst,ipk,mpk,ppk,0);
fprintf('err=%.4g\n',err)
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
    lr1=log(pa.r1o) + pa.r1e + pa.r1q;
    lm1=log(pa.m1o) + pa.m1e + pa.m1q;
    lrm1=lr1-lm1;
    if (lrm1>14) err=err+lrm1; end % limit tc1
    lr2=log(pa.r2o) + pa.r2e + pa.r2q;
    lm2=log(pa.m2o) + pa.m2e + pa.m2q;
    lrm2=lr2-lm2;
    if (lrm2>14) err=err+lrm2; end % limit tc1
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
    [~,ix]=max(D1(ii,k));
    mx=D1(ix,k);
    ph=unwrap(angle(Dh(ii,k)))/(2*pi);
    err=err+sum(ph>ph(1));                  % positive phase slope ?
    err=err+abs(mx-mpk(k))*2;               % peak
    err=err+abs(ph(ix)-ppk(k))*2;           % phase
    err=err+abs((ix-ipk(k)))*20;            % fp-map
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
k=k+1;pa.aco=pv(k);
k=k+1;pa.k1e=pv(k);
k=k+1;pa.r1e=pv(k);
k=k+1;pa.m1e=pv(k);
k=k+1;pa.k2e=pv(k);
k=k+1;pa.r2e=pv(k);
k=k+1;pa.m2e=pv(k);
k=k+1;pa.k3e=pv(k);
k=k+1;pa.r3e=pv(k);
k=k+1;pa.k4e=pv(k);
k=k+1;pa.ace=pv(k);
k=k+1;pa.k1q=pv(k);
k=k+1;pa.r1q=pv(k);
k=k+1;pa.m1q=pv(k);
k=k+1;pa.k2q=pv(k);
k=k+1;pa.r2q=pv(k);
k=k+1;pa.m2q=pv(k);
k=k+1;pa.k3q=pv(k);
k=k+1;pa.r3q=pv(k);
k=k+1;pa.k4q=pv(k);
k=k+1;pa.acq=pv(k);
end

%==============================================================

function [xx,Yb,Pd,Db,Dh,Zc,Gf]=fdmod23(pa,flst)
if (pa.m<1)
    xx=linspace(0,1,pa.n);
    Yb=xx;Pd=xx;Db=xx;Dh=xx;
    s=2i*pi*flst;
    Zc=pa.rco+pa.mco*s;
    ast=pa.ast;
    ama=pa.ama;
    gme=pa.gme;
    Zma=s*pa.mma+pa.rma+pa.kma./s;
    Zim=pa.rim+pa.kim./s;
    Zst=s*pa.mst+pa.rst+pa.kst./s;
    Gf=(ama*Zc)./(gme*(Zst+ast*Zc)+Zma.*(Zst+ast*Zc+Zim)./(gme*Zim));
    return;
end
L=pa.xl;
N=pa.n;
ast=pa.ast;
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
    srd=s*pa.rho*dx;
    [a,q,x,y,D]=mxfill(f,pa);    % initialize arrays
    p=mxsolve(a,q);
    [pd,yb,dd]=p_dif(p,y,f,D);
    [pd,yb,dd]=p_nan(pd,yb,dd);  % restrict plotting range
    Pd(:,k)=pd;
    Yb(:,k)=yb;
    Db(:,k)=dd(:,1);
    Dh(:,k)=dd(:,2);
    Vs=(p(nc,2)-p(nc,1))/srd;    % stapes linear velocity
    Zc(k)=p(nc,1)/(ast*Vs);      % cochlear impedance
end
xx=x/pa.xl;
Gf=Pd(1,:);
end

%--------------------------------------------------------------

function [a,q,x,y,D]=mxfill(f,pa)
% local parameter values
chsz=pa.chsz;
chsz=chsz*2/sum(chsz); % normalize channel sizes;
n=pa.n;
xl=pa.xl; rho=pa.rho; g=pa.gam;
% useful constructs
s=2i*pi*f;
dx=xl/(n-1);
srd=s*rho*dx;
% compute admittance for all x
x=transpose(linspace(0,xl,n));
[z1,z2,zh,za,ac,gh,bw]=imped(x,s,pa);
m=length(pa.chsz);
A=zeros(n,m);
Y=zeros(n,m,m);
a=zeros(3,m,m,n);
y=zeros(n,m,m);
% initialize apex
Zh=pa.khe/s+pa.rhe+pa.mhe*s;
z1(n)=Zh;z2(n)=1e-6;zh(n)=1e-6;za(n)=0;
% initialize partition admittance
zg=g*za;
zk=zeros(m,m);
if (m==1)
    D=1;                            % select PD
    B=2;                            % basal BC
elseif (m==2)
    D=[1 -1;0 0];                   % select PD
    B=[-1;1];                       % basal BC
elseif (m==3)
    D=[1 0 -1;0 [1 -1]; 0 0 0];     % select PD
    B=[-1;0;1];                     % basal BC
end
for k=1:n
    A(k,:)=ac(k)*chsz;
    if (m==1)
        hh=gh(k)*z2(k)/(z2(k)+zh(k));
        zk(1)=z1(k)+(gh(k)*zh(k)-zg(k))*hh;
    elseif (m==2)
        hh=gh(k)*z2(k)/(z2(k)+zh(k));
        zk(1,1)=z1(k)+(gh(k)*zh(k)-zg(k))*hh;
    elseif (m==3)
        zk(1,1:2)=[z1(k)  gh(k)*zh(k)-zg(k)]; 
        zk(2,1:2)=[z2(k)  -(zh(k)+z2(k))];
    end
    zk(m,:)=1; % conserve fluid volume
    Y(k,:,:)=zk\D;
    if (m==1)
        y(k,1,1)=1/zk(1,1);   % select Ybm 
        y(k,2,:)=hh*y(k,1,:); % select Yhb
    elseif (m==2)
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
% middle-ear boundary conditions
[alfx,qst]=midear(f,pa);
alfx=alfx*1e-7;
a(2,:,:,1)=A1*(1+alfx);
a(3,:,:,1)=A2*( -alfx);
% initialize q
q=zeros(m,n);
q(:,1)=qst*B;
end % return

function [z1,z2,z3,z4,ac,gh,bw]=imped(x,s,pa)
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
gh=pa.gpo*exp(pa.gpe*x+pa.gpq*q);
bw=pa.bwo*exp(pa.bwe*x+pa.bwq*q);
end

function [alfx,qst]=midear(f,pa)
% middle-ear impedances
s=2i*pi*f(:);
dx=pa.xl/(pa.n-1);
srd=s*pa.rho*dx;
ast=pa.ast;
ama=pa.ama;
gme=pa.gme;
if(pa.nmev==1)
    Zme=pa.mme*s+pa.rme+pa.kme./s;
    alfx=(Zme)./(srd*ast);
    qst =(ama/gme)/ast;
else
    Zma=s*pa.mma+pa.rma+pa.kma./s;
    Zim=pa.rim+pa.kim./s;
    Zst=s*pa.mst+pa.rst+pa.kst./s;
    gim=Zim./(Zma+pa.gme^2*Zim);
    alfx=(Zst+Zma.*gim)/(srd*ast)/1000; % <-- why ??
    qst =(ama*gme*gim)/50; % <-- why ??
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

function pa=modpar23(nch)
% chamber sizes
if (nch==1)
    exit(' single-chamber model not implemented.')
elseif (nch==3)
    pa=modpar24c3;
else
    pa=par_CEL16;
    %pa.gam=0.8;
end
pa.m=nch;             % number of fluid chambers
pa.xp=1;              % length of excitation pattern
end

%--------------------------------------------------------------

function pa=modpar24c3
pa.parlab='cel24c3';
pa.chsz=[0.9 0.1 1];           % three-chamber sizes
pa.gam = 1;                    % NDR multiplier
pa.n = 701;                    % number of points along BM
pa.xl = 3.5;                   % scala length
pa.yw = 0.1;                   % scala width
pa.zh = 0.1;                   % scala height
pa.rho = 1;                    % fluid density
pa.isv = [562 486 408 325 235 138 34]; % BM locations to save
pa.bwo = 0.05; pa.bwe = 0; pa.bwq = 0; % BM width
pa.gpo = 1.00; pa.gpe = 0; pa.gpq = 0; % partition gain (HB re BM)
pa.xtap=0; pa.xtex=6;
pa.hbt=0; pa.xp=1;
pa.khe=0.0001; pa.rhe=0.0001; pa.mhe=0.01;
% middle-ear parameters
pa.nmev=4; pa.kme=0.1; pa.rme=400; pa.mme=0.1;
pa.acp=0.5; pa.adi=0.2; pa.aed=0.5; pa.ama=0.5; pa.ast=0.03;
pa.cep=2700;  pa.rep=29.16; pa.mevgn=1.262; pa.stim=3;
pa.mcp=0.0002; pa.rcp=0.1; pa.kcp=2200; pa.rfz=0.01;
pa.mdi=0.005;  pa.rdi=100;  pa.kdi=7.7e+06;
pa.mrw=5e-3; pa.rrw=20; pa.krw=1.5e5; pa.arw=0.0625;
pa.mco=3; pa.rco=1e5;  pa.gme=0.5;
pa.kma=3e+06;  pa.rma=8000;  pa.mma=0.85;
pa.ked=7.04632e+06;  pa.red=200;  pa.med=0.02;
pa.kst=4.66022e+06;  pa.rst=8140.1;  pa.mst=0.89833;
pa.kim=3e+08;  pa.rim=28508.7; % nch=0 nme=4
% ---- parfit ----
pa.k1o=1.01219e+06;
pa.r1o=0.229798;
pa.m1o=2.58042e-05;
pa.k2o=302.994;
pa.r2o=30.4705;
pa.m2o=5.69822e-06;
pa.k3o=915257;
pa.r3o=0.00262073;
pa.k4o=237762;
pa.r4o=0;
pa.aco=0.0433444;
pa.gpo=0.493317;
pa.k1e=-2.4542;
pa.r1e=2.8997;
pa.m1e=-0.0202;
pa.k2e=-0.3602;
pa.r2e=-1.6061;
pa.m2e=-2.4805;
pa.k3e=-1.9948;
pa.r3e=2.3778;
pa.k4e=0.0000;
pa.r4e=0.0000;
pa.ace=0.0000;
pa.gpe=0.0000;
pa.k1q=0.222449;
pa.r1q=-1.304013;
pa.m1q=-0.056076;
pa.k2q=-0.031026;
pa.r2q=0.607868;
pa.m2q=0.195625;
pa.k3q=-0.000042;
pa.r3q=-0.003633;
pa.k4q=-0.801839;
pa.r4q=0.000000;
pa.acq=0.000000;
pa.gpq=0.000034;
pa.hbt=0.0000;
end

