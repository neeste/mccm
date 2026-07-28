; 2DOF.plt - mechanical elements
;
axlwt=0 : annlwt=0 : head=0
xllc=0 : yllc=0
xlen=10.5 : ylen=8
xmax=10.5 : ymax=8
mhal=c : mval=b : msgsiz=1
;
;----------------------- vertical mechanical elements ----------------
;
%define base ; <x1> <y1> <x2> <y2>
pltlwt=1 : pltype=rect : shade=14
($1) ($2)
($3) ($4)
plot
pltlwt=2 : pltype=line
($1) ($2)
($3) ($2)
plot
%%
%define mass ; <x1> <y1> <x2> <y2>
pltlwt=1 : pltype=rect : shade=0
($1) ($2)
($3) ($4)
plot
%%
; vertical spring
%define vspr ; <x1> <y1> <y2> <width> <height>
pltlwt=1 : pltype=line
_y1=($2+$3-$5)/2 : _y2=($2+$3+$5)/2
($1) ($2)
($1) (_y1)
rx=($4)/2 : ry=($4)/4 : _dy=($5-2*ry) : ndat=400 : tph=3.1415927*7
xdata=$$1+rx*sin(tph*(($$0-1)/(ndata+1)))
ydata=$$2-ry*cos(tph*(($$0-1)/(ndata+1)))+ry+_dy*(($$0-1)/(ndata+1))
($1) (_y1)
ndata=1 : xdata= : ydata=
($1) (_y2)
($1) ($3)
plot
%%
; vertical dashpot
%define vdsh ; <x1> <y1> <y2> <width> <height>
pltlwt=1 : pltype=line
_x1=($1-$4/2)    : _x2=($1+$4/2) 
_y1=($2+$3-$5)/2 : _y2=($2+$3+$5)/2
_x3=($1-$4/4)    : _x4=($1+$4/4) 
_y3=(_y1+_y2)/2
($1)  ($2)
($1)  (_y1)
plot
(_x1) (_y2)
(_x1) (_y1)
(_x2) (_y1)
(_x2) (_y2)
plot
(_x3) (_y3)
(_x4) (_y3)
;plot
($1)  (_y3)
($1)  ($3)
plot
%%
;
;---------------- 2DOF model --------------------
;
%base 3.0 5.0 6.3 5.3
%base 3.0 1.0 6.3 0.7
%mass 4.0 2.0 5.5 2.5
%mass 4.0 3.5 5.0 4.0
%vspr 4.3 1.0 2.0 0.2 0.5
%vspr 4.3 2.5 3.5 0.2 0.5
%vspr 4.3 4.0 5.0 0.2 0.5
%vdsh 4.7 1.0 2.0 0.3 0.2
%vdsh 4.7 2.5 3.5 0.3 0.2
%vdsh 4.7 4.0 5.0 0.3 0.2
; arrows
pltype=line : aronum=1 : arosiz=0.5
; M1
4.0 2.25
3.6 2.25
3.6 2.75
plot
; M2
4.0 3.75
3.6 3.75
3.6 4.25
plot
; Pd
5.4 1.4
5.4 2.0
plot
; Pa
5.4 3.1
5.4 2.5
plot
; labels
mval=3 : msgsiz=0.8
4.0  1.5  "k[1]"
4.0  3.0  "k[3]"
4.0  4.5  "k[2]"
5.1  1.5  "c[1]"
5.1  3.0  "c[3]"
5.1  4.5  "c[2]"
4.53 2.27 "m[1]"
4.51 3.77 "m[2]"
5.65 1.65 "P[d]"
5.65 2.85 "P[a]"
3.4  2.35 "^j^[1]"
3.4  3.85 "^j^[2]"
