; me_cmp.plt - compare middle-ear impedance & gain between measurements & model
;
head=0 : ticsiz=0.7 : ticdir=in : yhor=y
axlwt=0.7 : annlwt=1.1 : lablwt=1.1 : pltlwt=1.1 : labsiz=1.1 : sizfac=6
xmin=0.1 : xmax=10 : xcyc=1 : xper=100
xann=0.1 0.2 0.5 1 2 5 10
pltyp=line : clip=y
xdata=$1 : zdata=0
;---------------- zme ---------------------
;
;---------------- phase
;
newframe
xllc=1.5 : yllc=2.0
xlen=3 : ylen=1.5
xanskp=0
ymin=-0.25 : ymax=0.25 : ycyc=0 : yint=4 : yfmt=f.2 : yper=80
yann=-0.50 -0.25 0 0.25
xlabel=frequency (kHz)
ylabel=phase (cyc)
tlabel=
pltcol=1 : lintyp=0 : ydata=$2
include voss_phase.txt
plot
pltcol=2 : lintyp=0 : ydata=$2
include puria_phase.txt
plot
pltcol=4 : lintyp=0 : ydata=-$12
include yme_id.txt
plot
;---------------- magnitude
newframe
yllc=(yllc+ylen)
ylen=3
xanskp=-1
ymin=200 : ymax=5000 : ycyc=3 : yfmt=i : yper=100
yann=10^(3),10{8}
xlabel=
ylabel=|\ Z[me]\ | (mks ohm)
tlabel=
pltcol=1 : lintyp=0 : ydata=10^($2/10)
include voss_magnitude.txt
plot
pltcol=2 : lintyp=0 : ydata=10^($2/10)
include puria_magnitude.txt
plot
pltcol=4 : lintyp=0 : lindot=0 : ydata=1/$6
include yme_id.txt
plot
;---------------- key ---------------------
msgdat=n : mhal=0 : mval=0 : msgsiz=1 : msglwt=1.1 : msgang=0 : mhkey=3.3
1.1 2.8 "
|_0,pltcol=1| Voss
|_0,pltcol=2| Puria
|_0,pltcol=4| Lewis
"

;---------------- gain ---------------------
;
head=0 : ticsiz=0.7 : ticdir=in : yhor=y
axlwt=0.7 : annlwt=1.1 : lablwt=1.1 : pltlwt=1.1 : labsiz=1.1 : sizfac=6
xmin=0.1 : xmax=10 : xcyc=1 : xper=100
xann=0.1 0.2 0.5 1 2 5 10
pltyp=line : clip=y
xdata=$1 : zdata=0
;---------------- zme ---------------------
;
;---------------- phase
;
newframe
xllc=6 : yllc=2.0
xlen=3 : ylen=1.5
xanskp=0
ymin=-0.75 : ymax=0.25 : ycyc=0 : yint=4 : yfmt=f.25 : yper=80
yann=-0.50 -0.25 0 0.25
xlabel=frequency (kHz)
ylabel=phase (cyc)
tlabel=
pltcol=2 : lintyp=0 : ydata=$2/360
include me_trans_phase.txt
plot
pltcol=4 : lintyp=0 : ydata=$3
include gme_cmp.txt
plot
;---------------- key ---------------------
msgdat=n : mhal=0 : mval=0 : msgsiz=1.0 : msglwt=1.1 : msgang=0 : mhkey=3.3
0.2 1.9 "
"
;---------------- magnitude
newframe
yllc=(yllc+ylen)
ylen=3
xanskp=-1
ymin=0 : ymax=20 : yint=2.2 : yfmt=i : yper=100*2/3
yann=
xlabel=
ylabel=|\ P[v]\ /\ P[ec]\ |\  (dB)
tlabel=
pltcol=2 : lintyp=0 : ydata=$2
include me_trans_magnitude.txt
plot
pltcol=4 : lintyp=0 : lindot=0 : ydata=$2+1
include gme_cmp.txt
plot
;---------------- key ---------------------
msgdat=n : mhal=0 : mval=0 : msgsiz=1 : msglwt=1.1 : msgang=0 : mhkey=3.3
0.6 0.7 "
|_0,pltcol=2| Puria
|_0,pltcol=4| Lewis
"

