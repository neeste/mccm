; ldgrw.plt - loudness growth
;
head=0 : ticsiz=0.7 : ticdir=in : yhor=y
axlwt=0.7 : annlwt=1.1 : lablwt=1.1 : pltlwt=0.9 : labsiz=1.1 : sizfac=6
pltyp=line : clip=y : lintyp=0
xdata=$1 : zdata=0
;
newframe
xllc=1.5 : yllc=2.0
xlen=3.5 : ylen=3.5
xanskp=0
xmin=0    : xmax=100 : xint=5 : xper=100*5/6
ymin=0.01 : ymax=100 : ycyc=1 : yper=100
xlabel=stimulus level (dB SPL)
ylabel=loudness (sone)
tlabel=
;
; Fletcher & Munsen (1933, Table III)
lintyp=4 : pltlwt=0.5 : ydat=$2/975
include FM33.txt
plot
lintyp=0 : pltlwt=1 : ndat=1 : xdat=$1 : ydata=$2
;
pltcol=1 : ydata=$2
include ldngrwlv.txt
plot
;---------------- key ---------------------
msgdat=n : mhal=0 : mval=0 : msgsiz=1 : msglwt=1.1 : msgang=0 : mhkey=3.3
0.5 3.3 "
|_0,pltcol=1| tdm24
|_4,pltcol=0| Fletcher
"
;---------------- categorical loudness
newframe
xllc=6.2 : yllc=2.0
xanskp=0 : ycyc=0
xmin=0 : xmax=100 : xint=5 : xper=100*5/6
ymin=0 : ymax=100 : yint=5 : yper=100*5/6
xlabel=stimulus level (dB SPL)
ylabel=loudness (phon)
tlabel=
;
lintyp=4 : pltlwt=0.5 : pltcol=0.5 : ydata=$2
  0   0
100 100
plot
lintyp=0 : pltlwt=1   : pltcol=1   : ydata=$4
include ldngrwlv.txt
plot
;---------------- key ---------------------
msgdat=n : mhal=0 : mval=0 : msgsiz=1 : msglwt=1.1 : msgang=0 : mhkey=3.3
0.5 3.3 "
|_0,pltcol=1| tdm24
"
;
