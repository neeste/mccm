; moh24hbtc.plt - plot iso-displacement tuning curves
header=0 : clip=y : grid=0 : yhor=yes : ticdir=in
annsiz=1.6 : labsiz=1.6 : msgsiz=1.6 : sizfac=4
; part a
newframe
xllc=2   : xlen=4   : xcyc=3 : xmin=0.2 : xmax=20  : xans=-1 : xfmt=i
yllc=4   : ylen=3   : ycyc=0 : ymin=-20 : ymax=100 : yans=0 : yfmt=i : yint=6.2
yann=0 20 40 60 80 100
xlab=
ylab=threshold (dB SPL)
tlab=iso-displacement
zdata=0
pltype=lines : lintyp=0 : pltcol=1 : xdata=$1
%for 2:8
ydata=-$$$1
include moh24hbtc.txt
plot
%%
pltype=symbol : symbol=2 : symsiz=1 : pltcol=4 : ydata=$2
ydata=-$2
include moh24hbpf.txt
plot
symbol=11 : symsiz=0.5 : pltcol=1 : xdata=$4 : ydata=-$5
include moh24hbpf.txt
plot
;----------------
msgdat=n : mhal=0 : msgsiz=1.2 : mvkey=0.8
0.4 0.4 "
|2,pltcol=4| minimum audible pressure
"
;----------------
; part b
newframe
xllc=2 : xlen=4   : xcycle=3 : xmin=0.2 : xmax=20 : xanskp=0
yllc=2 : ylen=2   : ycycle=0 : ymin=-25  : ymax=5 : yint=6
xann=0.2 0.4 1 2 4 10 20
yann=-20 -10 0
xlab=frequency (kHz)
ylab=phase (cyc)
tlab=
pltype=lines only : lintyp=0 : pltcol=1 : xdata=$1
%for 9:15
ydata=$$$1
include moh24hbtc.txt
plot
%%
pltype=symbol : symbol=1 : symsiz=1 : pltcol=4 : ydata=$3
include moh24hbpf.txt
plot
symbol=11 : symsiz=0.5 : pltcol=1 : xdata=$4 : ydata=$6
include moh24hbpf.txt
plot
;----------------
msgdat=n : mhal=0 : msgsiz=1.2 : mvkey=0.8
2.0 (ylen-0.1) "
|1,pltcol=4| forward latency
"
