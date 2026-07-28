function [a,q,x,y,D]=mxfill(f,pa)
% set parameter values
m=pa.m; n=pa.n;
xl=pa.xl; rho=pa.rho; ast=pa.ast; ahe=pa.ahe; g=pa.gam;
aco=pa.aco; ace=pa.ace; bwo=pa.bwo; bwe=pa.bwe;
% useful constructs
s=2i*pi*f;
dx=xl/(n-1);
srd=s*rho*dx;
% compute admittance for all x
x=transpose(linspace(0,xl,n));
[z1,z2,zh,za,zoc]=imped(x,s,pa);
ac=aco*exp(ace*x);
bw=bwo*exp(bwe*x);
zg=zh-g*za;
y=zeros(n,2,m);
A=zeros(n,m);
Y=zeros(n,m,m);
a=zeros(3,m,m,n);
% initialize partition admittance
if (pa.tl)
   chsz=[1 1 1e9 1e9];
else
   chsz=[1 1 [1 1]*pa.mcsz];
end
D=[1 0 0 -1;0 0 0 0;0 0 0 0;0 0 0 0];
B=[1;-1;0;0];
H=[1 -1 0 0;-1 1 0 0;0 0 0 0;0 0 0 0];
for k=1:n
   if (pa.tl)
      zk=[z1(k)+zg(k) -zg(k) 0;
          -zh(k) z2(k)+zh(k) 0
              0      0      zoc(k)];
   else
      zk=[z1(k)          0  -zoc(k);
         -z2(k) z2(k)+zh(k)      0
             0        zg(k)  zoc(k)];
   end
   yp=zk\[1;0;0];
   zk(:,m)=0;
   zk(m,:)=1;
   A(k,:)=ac(k)*chsz;
   Y(k,:,:)=zk\D;
   y(k,1,1)=yp(1);
   y(k,2,1)=yp(1)-yp(2);
   if (pa.tl==0)
      y(k,1,1)=yp(1);
      y(k,2,1)=yp(2);
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
