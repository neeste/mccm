; me_circuit.plt - ME network with three KRM branches + Rco + CR + ISJ
;
magn=1.2  : wylen=8 : ss=0.5 : lw=0.8
axlwt=0 : annlwt=0 : head=0
xllc=0 : yllc=0
xlen=10.5 : ylen=8
xmax=10.5 : ymax=8
;
;----------------------- circuit elements ----------------
;
; terminals
%define term ; <x1> <y1> <x2> <y2>
pltlwt=lw : pltype=symbols : symsiz=ss : symbol=1
($1) ($2)
($3) ($2)
plot
($1) ($4)
($3) ($4)
plot
%%
; wires
%define wire ; <x1> <y1> <x2> <y2>
pltlwt=lw : pltype=line
($1) ($2)
($3) ($4)
plot
%%
;
;----------------------- vertical circuit elements ----------------
;
; resistor - vertical
%define vres ; <x1> <y1> <y2> <width> <height>
pltlwt=lw : pltype=line
_dy=($5)/7 : _dx=($4)/2
_x1=$1+_dx : _x2=$1-_dx
_y1=($2+$3-$5-_dy)/2 : _y2=($2+$3+$5-_dy)/2
($1)      ($2)
($1)      (_y1)
(_x2) (_y1+_dy*1)
(_x1) (_y1+_dy*2)
(_x2) (_y1+_dy*3)
(_x1) (_y1+_dy*4)
(_x2) (_y1+_dy*5)
(_x1) (_y1+_dy*6)
($1)      (_y2)
($1)      ($3)
plot
%%
; capacitor - vertical
%define vcap ; <x1> <y1> <y2> <width> <height>
pltlwt=lw : pltype=line
_x1=($1-$4/2)    : _x2=($1+$4/2) 
_y1=($2+$3-$5)/2 : _y2=($2+$3+$5)/2
($1)  ($2)
($1)  (_y1)
plot
(_x1) (_y1)
(_x2) (_y1)
plot
(_x1) (_y2)
(_x2) (_y2)
plot
($1)  (_y2)
($1)  ($3)
plot
%%
%%; inductor - vertical
%define vind ; <x1> <y1> <y2> <width> <height>
pltlwt=lw : pltype=line
_y1=($2+$3-$5)/2 : _y2=($2+$3+$5)/2
($1) ($2)
($1) (_y1)
rx=($4)/2 : ry=($4)/5 : _dy=($5-2*ry) : ndat=400 : tph=3.1415927*7
xdata=$$1+rx*sin(tph*(($$0-1)/(ndata+1)))
ydata=$$2-ry*cos(tph*(($$0-1)/(ndata+1)))+ry+_dy*(($$0-1)/(ndata+1))
($1) (_y1)
ndata=1 : xdata= : ydata=
($1) (_y2)
($1) ($3)
plot
%%
%%; transfmormer - vertical
%define vtrn ; <x1> <x2> <y1> <y2> <width> <height>
pltlwt=lw : pltype=line
rx=($5)/2 : ry=($5)/5 : _dy=($6-2*ry) : tph=6.2832*3.5
_x1=($1) : _y1=($3+$4-$6)/2 : _y2=($3+$4+$6)/2 : _dx=($2-$1-$5*5/4)/2
(_x1) ($3)
(_x1) (_y1)
ndata=400
xdata=$$1+rx*sin(tph*(($$0-1)/(ndata+1)))
ydata=$$2-ry*cos(tph*(($$0-1)/(ndata+1)))+ry+_dy*(($$0-1)/(ndata+1))
(_x1+_dx) (_y1)
ndata=1 : xdata= : ydata=
(_x1) (_y2)
(_x1) ($4)
plot
_x1=($2) : _y1=($3+$4-$6)/2 : _y2=($3+$4+$6)/2
(_x1) ($3)
(_x1) (_y1)
ndata=400
xdata=$$1-rx*sin(tph*(($$0-1)/(ndata+1)))
ydata=$$2-ry*cos(tph*(($$0-1)/(ndata+1)))+ry+_dy*(($$0-1)/(ndata+1))
(_x1-_dx) (_y1)
ndata=1 : xdata= : ydata=
(_x1) (_y2)
(_x1) ($4)
plot
%%
; transmission line - vertical
%define htln ; <x1> <y1> <x2> <y2>
pltlwt=2 : pltype=line
($1) ($2)
($1) ($4)
plot
($3) ($2)
($3) ($4)
plot
pltlwt=lw : pltype=symbols : symsiz=ss
($1) ($2)
($1) ($4)
plot
($3) ($2)
($3) ($4)
plot
%%
; source - vertical
%define vsrc ; <x1> <y1> <y2> <width> <height>
pltlwt=lw : pltype=line
($1) ($2)
($1) ($3)
plot
pltlwt=lw : pltype=symbols : symsiz=2 : symbol=1
($1) ($2+$3)/2
plot
msgsiz=1.2
($1) ($2+$3-0.4)/2 "~"
msgsiz=1
%%
;
;----------------------- horizontal circuit elements ----------------
;
; resistor - horizontal
%define hres ; <x1> <x2> <y1> <width> <height>
pltlwt=lw : pltype=line
_dx=($4)/7 : _dy=($5)/2
_x1=($1+$2-$4-_dx)/2 : _x2=($1+$2+$4-_dx)/2
_y1=($3-_dy) : _y2=($3+_dy)
($1)        ($3) 
(_x1)       ($3)
(_x1)       ($3)    
(_x1+_dx*1) (_y2) 
(_x1+_dx*2) (_y1) 
(_x1+_dx*3) (_y2) 
(_x1+_dx*4) (_y1) 
(_x1+_dx*5) (_y2) 
(_x1+_dx*6) (_y1) 
(_x2)       ($3)
($2)        ($3) 
plot
%%
; capacitor - horizontal
%define hcap ; <x1> <x2> <y1> <width> <height>
pltlwt=lw : pltype=line
_x1=($1+$2-$4)/2 : _x2=($1+$2+$4)/2
_y1=($3-$5/2) : _y2=($3+$5/2)
($1)  ($3) 
(_x1) ($3)
plot
(_x1) (_y1)
(_x1) (_y2)
plot
(_x2) (_y1)
(_x2) (_y2)
plot
(_x2) ($3)
($2)  ($3) 
plot
%%
; inductor - horizontal
%define hind ; <x1> <x2> <y1> <width> <height>
pltlwt=lw : pltype=line
_x1=($1+$2-$4)/2 : _x2=($1+$2+$4)/2 
($2)  ($3) 
(_x2) ($3)
rx=($5)/5 : ry=($5)/2 : _dx=($4-2*rx) : ndat=400 : tph=6.2832*3.5
xdata=rx*cos(tph*(($$0-1)/(ndata+1)))+$$1-rx-_dx*(($$0-1)/(ndata+1))
ydata=ry*sin(tph*(($$0-1)/(ndata+1)))+$$2
(_x2) ($3)
ndata=1 : xdata= : ydata=
(_x1) ($3)
($1)  ($3) 
plot
%%
; transmission line - horizontal
%define htln ; <x1> <y1> <x2> <y2>
pltlwt=2 : pltype=line
($1) ($2)
($3) ($2)
plot
($1) ($4)
($3) ($4)
plot
pltlwt=lw : pltype=symbols : symsiz=ss : symbol=1
($1) ($2)
($3) ($2)
plot
($1) ($4)
($3) ($4)
plot
%%
;=================================================================================
;
;
;---------------- three KRM branches --------------------
;
o=3
%hcap 6.1 6.6 o+1.5 0.05 0.2
%hres 6.6 7.1 o+1.5 0.30 0.2
%hind 7.1 7.6 o+1.5 0.30 0.2
%wire 7.6 o+1.5 7.8 o+1.5
%vres 7.8 o+0.0 o+1.5 0.2 0.30
%hcap 4.9 5.6 o+0.5 0.05 0.2
%hres 4.9 5.6 o+0.0 0.30 0.2
%hind 4.9 5.6 o-0.5 0.30 0.2
%vcap 6.0 o+0.75 o+1.2 0.2 0.05
%vres 6.0 o+0.3 o+0.75 0.2 0.30
%wire 6.0 o+0.0 6.0 o+0.3
%wire 6.0 o+1.2 6.0 o+1.5
%vcap 4.1 o+1.0 o+1.5 0.2 0.05
%vres 4.1 o+0.6 o+1.1 0.2 0.30
%vind 4.1 o+0.0 o+0.6 0.2 0.30
%hcap 4.4 4.9 o+1.5 0.05 0.2
%hres 4.9 5.4 o+1.5 0.30 0.2
%hind 5.4 5.9 o+1.5 0.30 0.2
%wire 5.9 o+1.5 6.1 o+1.5
%wire 3.4 o+1.5 3.4 o+1.0
%wire 3.4 o+1.5 4.4 o+1.5
%wire 3.4 o+0.5 3.4 o+0.0
%wire 3.4 o+0.0 4.9 o+0.0
%wire 5.6 o+0.0 7.8 o+0.0
%wire 4.9 o-0.5 4.9 o+0.5
%wire 5.6 o-0.5 5.6 o+0.5
%htln 1.4 o+0.5 3.4 o+0.5
%htln 1.4 o+1.0 3.4 o+1.0
;
; non-uniform transmission line
xdata=$1+$0*2/6
ydata=$2
ndata=5
symbol=1 : symsiz=0.3
1.4 o+0.5
plot
1.4 o+1.0
plot
xdata=
ydata=
ndata=1
; ------------------------------
mhal=l : mval=c : msgsiz=0.7 : msgang=0
1.0 o+0.85 "Z[ec]"
3.4 o+0.85 "Z[me]"
mhal=r : mval=b : msgsiz=0.6
6.50 o+1.70 "K[1]"
7.00 o+1.70 "R[st]"
7.50 o+1.70 "M[1]"
7.60 o+0.75 "R[co]"
4.55 o+1.20 "K[2]"
4.55 o+0.80 "R[2]"
4.55 o+0.30 "M[2]"
4.80 o+1.70 "K[3]"
5.30 o+1.70 "R[3]"
5.80 o+1.70 "M[3]"
5.40 o+0.68 "K[4]"
5.40 o+0.18 "R[4]"
5.40 o-0.32 "M[4]"
6.45 o+0.90 "K[5]"
6.45 o+0.45 "R[5]"
mhal=c : mval=b : msgsiz=0.6
2.4 o+1.1 "ear canal"
mhal=c : mval=b : msgsiz=0.6 : msgang=90
;0.4 o+0.8 "stapes\ndisarticulation"
;
