; moh24tc.plt - plot iso-displacement tuning curves
header=0 : clip=y : grid=0 : yhor=yes : ticdir=in
annsiz=1.6 : labsiz=1.6 : msgsiz=1.6 : sizfac=4 : solid=0
; part a
newframe
xllc=2 : xlen=4 : xcyc=0 : xmin=0 : xmax=1  : xans=-1 : xfmt=f.1
yllc=4 : ylen=3 : ycyc=0 : ymin=0 : ymax=100 : yans=0 : yper=100*5/6
yann=0 20 40 60 80 100
xlab=
ylab=magnitude (dB)
tlab=excitation patterns
zdata=0
pltype=lines : lintyp=0 : pltcol=1 : xdata=$1
%for 2:8
ydata=100+$$$1
include moh24ep.txt
plot
%%
pltype=symbol : symbol=1 : symsiz=1 : pltcol=4 : ydata=$2
ydata=100+$2
include moh24px.txt
plot
;symbol=11 : symsiz=0.5 : pltcol=1 : xdata=$4 : ydata=100+$5
;include moh24px.txt
;plot
; part b
newframe
xllc=2 : xlen=4 : xcycle=0 : xmin=0 : xmax=1 : xanskp=0
yllc=2 : ylen=2 : ycycle=0 : ymin=-25  : ymax=5 : yint=6 : yper=100
yann=-20 -10 0
xlab=x / L
ylab=phase (cyc)
tlab=
pltype=lines only : lintyp=0 : pltcol=1 : xdata=$1
%for 9:15
ydata=$$$1
include moh24ep.txt
plot
%%
pltype=symbol : symbol=1 : symsiz=1 : pltcol=4 : ydata=$3
include moh24px.txt
plot
;symbol=11 : symsiz=0.5 : pltcol=1 : xdata=$4 : ydata=$6
;include moh24px.txt
;plot
; -------- labels
