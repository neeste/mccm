% fit parameters to reference set

function fitref
f=1000*2.^(-1:3);
% fetch model parameters
pa=modpar16;pa.m=4;
[x0,cf0,d0,~]=fetch_set(f,pa);
pa=modpar17;pa.m=4;
[x1,cf1,d1,pd]=fetch_set(f,pa);
xlim=[0 pa.xl];
plot_sets(x0,cf0,d0,x1,cf1,d1,pd,xlim);
% search for optimal parameter values
op=optimset('MaxFunEvals',2000);
for k=1:10
   pv=get_values(pa);
   pv9=pv(1:9);
   pv(pv9<1e-9)=1e-9;
   err1=refdev(pv,f,pa,d0,cf0);
   if (k==1)err0=err1;end
   pv=fminsearch(@(pv) refdev(pv,f,pa,d0,cf0),pv,op);
   pa=set_values(pa,pv);
   pa.gmo=1;
   err2=refdev(pv,f,pa,d0,cf0);
   fprintf('fitpar(%d): err=%6.3f,%6.3f,%6.3f\n',k,err0,err1,err2);
   [x1,cf1,d1,pd]=fetch_set(f,pa);
   plot_sets(x0,cf0,d0,x1,cf1,d1,pd,xlim);
   wrtpar(pa,'fitref.txt');
end
return

%------------------------------------------------------

function plot_sets(x0,cf0,d0,x1,cf1,d1,pd,xlim)
% plot data
[~,nf,ng]=size(d0);
lab='HB displacement (nm) at 0 dB SPL';
for j=1:ng
   mdxfrm(j,xlim,[-120 -20],[-5 20],lab);
   for k=1:nf
      mdxplt(j,x0,cf0(:,j),d0(:,k,j),x1,cf1(:,j),d1(:,k,j));
   end
end
%
latency(ng+1,cf0,d0,cf1,d1)
%
lab='BM presure difference';
mpfrm(ng+2,[xlim  -70  30],[xlim -9.0 1.0],lab);
for k=1:nf
   mpplt(ng+2,x1,pd(:,k));
end
drawnow;
return

function latency(fig,o0,d0,o1,d1)
[~,nf,ng]=size(d0);
td0=zeros(nf,ng);
td1=zeros(nf,ng);
fd0=zeros(nf,ng);
fd1=zeros(nf,ng);
for j=1:ng
   cf0=2.^o0(:,j);
   cf1=2.^o1(:,j);
   for k=1:nf
      mg0=20*log10(max(1e-9,abs(d0(:,k,j))));
      mg1=20*log10(max(1e-9,abs(d1(:,k,j))));
      ph0=unwrap(angle(d0(:,k,j)))/(2*pi);
      ph1=unwrap(angle(d1(:,k,j)))/(2*pi);
      gd0=cdif(ph0)./cdif(cf0);
      gd1=cdif(ph1)./cdif(cf1);
      [~,xx0]=max(mg0);
      [~,xx1]=max(mg1);
      td0(k,j)=gd0(xx0);
      td1(k,j)=gd1(xx1);
      fd0(k,j)=cf0(xx0);
      fd1(k,j)=cf1(xx1);
   end
end
% plot
gdlab='latency (ms)';
x_lab='frequency (kHz)';
figure(fig);clf
loglog(fd0,td0(:,1),'b',fd1,td1(:,1),'go-');
hold on % workaround to fix legend ???
loglog(fd0,td0,'b',fd1,td1,'go-');
hold off
axis([0.25 16 1 20]);
xlabel(x_lab);
ylabel(gdlab);
legend('ref','new');
return

function mxfrm(n,lim1,lim2,lim3,lab1,lab2,lab3,lab4)
figure(n);clf
subplot(211);hold on;axis([lim1 lim2]);
title(lab1);
ylabel(lab2);
subplot(212);hold on;axis([lim1 lim3]);
ylabel(lab3);
xlabel(lab4);
return

function mdxfrm(n,lim1,lim2,lim3,tlab)
mglab='magnitude (dB)';
gdlab='delay (ms)';
x_lab='distance from stapes (cm)';
mxfrm(n,lim1,lim2,lim3,tlab,mglab,gdlab,x_lab);
return

function mdxplt(n,x1,o1,y1,x2,o2,y2)
mg1=20*log10(max(1e-9,abs(y1)));
mg2=20*log10(max(1e-9,abs(y2)));
ph1=unwrap(angle(y1))/(2*pi);
ph2=unwrap(angle(y2))/(2*pi);
gd1=cdif(ph1)./cdif(2.^o1);
gd2=cdif(ph2)./cdif(2.^o2);
figure(n);
subplot(211);plot(x1,mg1,'b',x2,mg2,'g');
subplot(212);plot(x1,gd1,'b',x2,gd2,'g');
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

function pv=get_values(pa)
pv=zeros(20,1);
pv( 1)=pa.k1o;
pv( 2)=pa.k2o;
pv( 3)=pa.k3o;
pv( 4)=pa.k4o;
pv( 5)=pa.r1o;
pv( 6)=pa.r2o;
pv( 7)=pa.r3o;
pv( 8)=pa.m1o;
pv( 9)=pa.m2o;
pv(10)=pa.k1e;
pv(11)=pa.k2e;
pv(12)=pa.k3e;
pv(13)=pa.k4e;
pv(14)=pa.r1e;
pv(15)=pa.r2e;
pv(16)=pa.r3e;
pv(17)=pa.m1e;
pv(18)=pa.m2e;
pv(19)=pa.k5o;
pv(20)=pa.k5e;
pv(21)=pa.r5o;
pv(22)=pa.r5e;
%pv(23)=pa.r4o;
%pv(24)=pa.r4e;
return

function pa=set_values(pa,pv)
pa.k1o=pv( 1);
pa.k2o=pv( 2);
pa.k3o=pv( 3);
pa.k4o=pv( 4);
pa.r1o=pv( 5);
pa.r2o=pv( 6);
pa.r3o=pv( 7);
pa.m1o=pv( 8);
pa.m2o=pv( 9);
pa.k1e=pv(10);
pa.k2e=pv(11);
pa.k3e=pv(12);
pa.k4e=pv(13);
pa.r1e=pv(14);
pa.r2e=pv(15);
pa.r3e=pv(16);
pa.m1e=pv(17);
pa.m2e=pv(18);
pa.k5o=pv(19);
pa.k5e=pv(20);
pa.r5o=pv(21);
pa.r5e=pv(22);
%pa.r4o=pv(23);
%pa.r4e=pv(24);
return

%==============================================================

function [x,cf,dhb,pbm]=fetch_set(f,pa)
n=pa.n;
nf=length(f);
dhb=zeros(n,nf);
pbm=zeros(n,nf);
%
for k=1:nf
   [a,q,x,y,D]=mxfill(f(k),pa);
   p=mxsolve(a,q);
   [pd,~,dd,~]=p_dif(p,y,f(k),D);
   dhb(:,k)=dd(:,2);
   pbm(:,k)=pd(:,1);
end
cf=fpmap(x,f,dhb);
%
return

function [cf,cp,c1]=fpmap(x,f,d)
nf=size(d,2);
mg=20*log10(max(1e-9,abs(d)));
cp=zeros(1,nf);
for k=1:nf
   [~,i]=max(mg(:,k));
   cp(k)=x(i);
end
fo=log2(f/1000);
c1=polyfit(fo(:),cp(:),1);
cf=(x-c1(2))/c1(1);
return

function d=cdif(x)
d=zeros(size(x));
n=length(x);
d(1)=x(2)-x(1);
d(n)=x(n)-x(n-1);
d(2:(n-1))=x(3:n)-x(1:(n-2));
return

%==============================================================

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

%==============================================================

function err=refdev(pv,f,pa,d0,of0)
[~,nf]=size(d0);
if(min(pv(1:9))<eps)err=1e9;return;end
pa=set_values(pa,pv);
if(pa.r3e<-4)err=1e9;return;end
[~,of1,d1,pd]=fetch_set(f,pa);
mxmg=zeros(1,nf);
mxgd=zeros(1,nf);
mngd=zeros(1,nf);
mxgd0=zeros(1,nf);
mxgd1=zeros(1,nf);
err=0;
cnt=0;
%
cf0=2.^of0;
cf1=2.^of1;
for k=1:nf
   y0=d0(:,k);
   y1=d1(:,k);
   mg0=20*log10(max(1e-9,abs(y0)));
   mg1=20*log10(max(1e-9,abs(y1)));
   ph0=unwrap(angle(y0))/(2*pi);
   ph1=unwrap(angle(y1))/(2*pi);
   gd0=cdif(ph0)./cdif(cf0);
   gd1=cdif(ph1)./cdif(cf1);
   [~,ix0]=max(mg0);
   [~,ix1]=max(mg1);
    mxmg(k)=mg1(ix1);
    mxgd(k)=gd1(ix1);
    mngd(k)=gd1(ix1);
    mxgd0(k)=gd0(ix0);
    mxgd1(k)=gd1(ix1);
    if(mxgd1(k)<eps)err=1e9;return;end
    err1=abs(log2(mxgd1(k)/mxgd0(k)))*80; % delay error
    d0k=d0(:,k);
    d1k=d1(:,k);
    ii=~(isnan(d0k)|isnan(d1k));
    a0k=sqrt(mean(abs(d0k(ii)).^2));
    a1k=sqrt(mean(abs(d1k(ii)).^2));
    ddk=d1k(ii)/a1k-d0k(ii)/a0k;
    err2=mean(abs(ddk).^2)*8;      % displacement error
    err3=mean((of0-of1).^2)*160;   % fpmap error
    ii=(ix0-1):(ix0+1);
    md=mg1(ii)-mg0(ii);
    err4=mean(abs(md));            % maxmag error
    php=unwrap(angle(pd(:,k)))/(2*pi);
    mxph=max(0,max(php));
    err5=mxph*20;                  % maxphp > 0
    err=err+err1+err2+err3+err4+err5;
    cnt=cnt+1;
end
dgd=diff(log(mxgd1));              % slope of log(gd) vs log(f)
err1=sum(dgd(dgd>0));              % penalize positive slope ???
err2=sum(abs(diff(dgd)))*1;        % penalize curvature
err=err+err1+err2;
err=err/cnt;
%
%err=err+log10(2*pa.r2o/pa.r3o)*80;  % penalize large R2/R3 ratio
%
[~,k4]=min(abs(f-4000));
[~,k8]=min(abs(f-8000));
err=err+abs(mxmg(k4)+35);         % target lo-lev 4kHz HB disp.
err=err+abs(mxmg(k8)+35);         % target lo-lev 8kHz HB disp.
fprintf('refdev: err=%6.3f\n',err);
return
