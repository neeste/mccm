; ldsum.plt - loudness summation (oct)
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
xmin=0 : xmax=100 : xint=5 : xper=100*5/6
ymin=0 : ymax=100 : yint=5 : yper=100*5/6
xlabel=stimulus level (dB SPL)
ylabel=loudness (phon)
tlabel=loudness growth
;
; CLS data - Rasetshwane et al. 2015
lintyp=1 : pltlwt=0.7 : pltype=both : symbol=11 : symsiz=0.3
pltcol=0 : ydata=$3
include ldnclsnh.txt
plot
;
pltyp=lin : lintyp=0 : pltlwt=0.9 : pltype=lines
pltcol=1 : ydata=$2
include ldnsumlv.txt
plot
pltcol=2 : ydata=$4
include ldnsumlv.txt
plot
pltcol=4 : ydata=$7
include ldnsumlv.txt
plot
;---------------- key ---------------------
msgdat=n : mhal=0 : mval=0 : msgsiz=1 : msglwt=1.1 : msgang=0 : mhkey=3.3
pltcol=0 : mvkey=0.3
0.2 3.3 "
|11||_1,pltcol=0,pltlwt=0.7||11|  Rasetshwane
|_0,pltcol=1,pltlwt=0.9|   0
|_0,pltcol=2,pltlwt=0.9|  1/4 octave
|_0,pltcol=4,pltlwt=0.9|   1  
"
2.5 1.0 "1 kHz"
;
newframe
xllc=6
xanskp=0
xmin=0.05 : xmax=1 : xcyc=1 : xper=100 : xfmt=f.2
ymin=50 : ymax=65  : yint=3.5   : yper=100
xann=0.1 0.2 0.5 1
xlabel=stimulus bandwidth (kHz)
ylabel=equivalent level (dB SPL)
tlabel=loudness summation
xdata=$1/1000 : zdata=0
;
; critical-band reference point
pltyp=symb : symb=7 : symsiz=1 : pltlwt=0.5
pltcol=8 : ydata=$2
133 60
plot
;
; Loudness-summation data - Leibold et al. 2007
lintyp=1 : pltlwt=0.7 : pltype=both : symbol=11 : symsiz=0.3
pltcol=0 : ydata=$2
  46 56.1
  92 55.3
 231 57.8
 465 60.3
 956 63.5
2119 66.8
plot
;
xdata=$1
pltype=lines : lintyp=0 : pltlwt=0.9
pltcol=0 : select= : ydata=$5
include ldnsumbw.txt
plot
pltype=symb : symbol=11 : symsiz=0.8
pltcol=2 : select=($1>0.16)&($1<0.26)
include ldnsumbw.txt
plot
pltcol=4 : select=$1>0.7
include ldnsumbw.txt
plot
;---------------- key ---------------------
msgdat=n : mhal=0 : mval=0 : msgsiz=1 : msglwt=1.1 : msgang=0 : mhkey=3.3
pltcol=0 : mvkey=0.3
0.2 3.3 "
|11||_1,pltcol=0,pltlwt=0.7||11|  Leibold 2007
|_0,pltcol=1,pltlwt=0.9|  model
"
