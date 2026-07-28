% two-chamber cochlear model

function admit
L=3.5;
mpfrm(1,[0 L -130 -30],[0 L -0.5 0.6],'BM admittance');
mpfrm(2,[0 L -60   40],[0 L -0.3 0.7],'HB transfer');
% loop over 10 frequencies
for fi=0:9
   f=400*(sqrt(2)^fi);
   % initialize arrays
   [x,y,h]=mxfill(f);
   % plot data
   mpplt(1,x,y(:,1));
   mpplt(2,x,h(:,1));
   drawnow;
end
return

function [x,y,h]=mxfill(f)
% set parameter values
n=701; xl=3.5;g=0;
% compute admittance for all x
x=transpose(linspace(0,xl,n));
s=2i*pi*f;
[z1,z2,zc,za]=imped(x,s);
w=[1 1];
y=zeros(n,2);
for k=1:n
    zack=zc(k)-g*za(k);
    zk=[z1(k)+zack -zack;-zc(k) z2(k)+zc(k)];
    wk=diag(w);
    y(k,:)=[1 0]*(zk\wk);
end
h=z2./(z2+zc); % HB
return

function [z1,z2,zc,za]=imped(x,s)
% set parameter values
k1o=1e9; r1o=800; m1o=0.15; k1e=-1.2; r1e=-0.2; m1e=0.8; w1=1;
k2o=4e6; r2o=4;  m2o=0.003; k2e=-1.2; r2e=-0.2; m2e=0.8; w2=1;
kco=1e7; rco=16;            kce=-1.2; rce=-0.2;
kao=3.6e7; rao=48;          kae=-1.2; rae=-0.2;
% compute impedance for all x
z1=k1o*exp(k1e*x)/s+r1o*exp(r1e*x)+m1o*exp(m1e*x)*s;
z2=k2o*exp(k2e*x)/s+r2o*exp(r2e*x)+m2o*exp(m2e*x)*s;
zc=kco*exp(kce*x)/s+rco*exp(rce*x);
za=kao*exp(kae*x)/s+rao*exp(rae*x);
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
mg=20*log10(max(1e-9,abs(y)));
ph=unwrap(angle(y))/(2*pi);
subplot(211);plot(x,mg);
subplot(212);plot(x,ph);
return
