; ldnnen.plt - loudness summation (Hz) across number of ensembles
;
head=0 : ticsiz=0.7 : ticdir=in : yhor=y
axlwt=0.7 : annlwt=1.1 : lablwt=1.1 : pltlwt=0.9 : labsiz=1.1 : sizfac=6
pltyp=line : clip=y : lintyp=0
xdata=$1/1000 : zdata=0
;
newframe
xllc=1.5 : yllc=2.0
xlen=3.5 : ylen=3.5
xanskp=0
xmin=0.05 : xmax=1 : xcyc=1 : xper=100
ymin=50   : ymax=65  : yint=3 : yper=100
xann=0.1 0.2 0.5 1
xlabel=stimulus bandwidth (kHz)
ylabel=equivalent loudness (dB SPL)
tlabel=ensemble number
;
; critical-band reference point
pltyp=symb : symb=7 : symsiz=1
pltcol=0 : ydata=$2
 133 60
plot
;
; Loudness Summation data - Leibold et al. 2007
lintyp=1 : pltlwt=0.5 : pltyp=both : symb=11 : symsiz=0.3
pltcol=0 : ydata=$2
  46 56.1
  92 55.3
 231 57.8
 465 60.3
 956 63.5
2119 66.8
plot
mvkey=0.3
0.2 0.3 "|11||_1||11| Leibold et al. (2007)"
pltyp=line : lintyp=0 : pltlwt=1 : a=16.4
;
xdata=$1
pltcol=1 : ydata=$2
include ldnnen.txt
plot
pltcol=2 : ydata=$3
include ldnnen.txt
plot
pltcol=4 : ydata=$4
include ldnnen.txt
plot
pltcol=3 : ydata=$5
include ldnnen.txt
plot
;---------------- key ---------------------
msgdat=n : mhal=0 : mval=0 : msgsiz=1 : msglwt=1.1 : msgang=0 : mhkey=3.3
0.3 3.3 "
|_0,pltcol=1|  35
|_0,pltcol=2|  70
|_0,pltcol=4| 140
"
