; moh24hbtc.plt - plot iso-displacement tuning curves
header=0 : clip=y : grid=0 : yhor=yes : ticdir=in
annsiz=1.6 : labsiz=1.6 : sizfac=4
; part a
newframe
xllc=2   : xlen=5   : xcyc=3 : xmin=0.1 : xmax=20  : xans=0 : xfmt=i
yllc=2   : ylen=5   : ycyc=0 : ymin=0   : ymax=100 : yans=0 : yfmt=i : yint=10 : yper=100
yann=0 20 40 60 80 100
xann=0.2 0.4 1 2 4 10 20
xlab=frequency (kHz)
ylab=threshold (dB SPL)
tlab=iso-displacement
zdata=0
pltype=lines : lintyp=0 : xdata=$1 : ydata=-$7
pltcol=1 : pltlwt=1.4
include moh24hbtc.txt
plot
pltcol=4 : pltlst=0.7
include moh24bmtc.txt
plot
;----------------
msgdat=n : mhal=0 : msgsiz=1.4 : mhkey=4
0.5 0.9 "
|_0,pltcol=1,pltlwt=1.4| hair bundle
|_0,pltcol=4,pltlwt=0.7| basilar membrane
"

