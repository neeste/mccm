; ldnsum.plt - loudness summation (oct): comparison with repetiions
;
head=0 : ticsiz=0.7 : ticdir=in : yhor=y
axlwt=0.7 : annlwt=1.1 : lablwt=1.1 : pltlwt=0.9 : labsiz=1.1 : sizfac=6
pltyp=line : clip=y : lintyp=0
xdata=$1 : zdata=0
;
newframe
xllc=1.5 : yllc=2.0
xlen=4.5 : ylen=4.5
;
newframe
xmin=0.05 : xmax=1 : xcyc=1   : xper=100 : xfmt=f.2
ymin=45 : ymax=65  : yint=4.5 : yper=100
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
; five-tone loudness-summation data - Leibold et al. 2007
lintyp=4 : pltlwt=0.7 : pltype=both : symbol=11 : symsiz=0.3
pltcol=1 : xdata=$1/1000
  46 56.1
  92 55.3
 231 57.8
 465 60.3
 956 63.5
2119 66.8
plot
; flat-noise loudness-summation data - CLS 20241126
pltcol=2 : xdata=$1
          0.005             60             60
         0.1735         52.869         55.604
        0.70711         55.583         58.152
plot
;
xdata=$2 : ydata=$5
include qolr_tone/ldnsumbw.txt
y1=60-ydat(1)
data
include qolr_flat/ldnsumbw.txt
y2=60-ydat(1)
data
include qolr_lono/ldnsumbw.txt
y3=60-ydat(1)
data
include qolr_amod/ldnsumbw.txt
y4=60-ydat(1)
data
;
xdata=$1
pltype=lines : lintyp=0 : pltlwt=0.9
;----------------
pltcol=1 : select= : ydata=$5
include qolr_tone/ldnsumbw1.txt
plot
include qolr_tone/ldnsumbw2.txt
plot
include qolr_tone/ldnsumbw3.txt
plot
include qolr_tone/ldnsumbw4.txt
plot
;----------------
pltcol=2 : select= : ydata=$5
include qolr_flat/ldnsumbw1.txt
plot
include qolr_flat/ldnsumbw2.txt
plot
include qolr_flat/ldnsumbw3.txt
plot
include qolr_flat/ldnsumbw4.txt
plot
;----------------
pltcol=4 : select= : ydata=$5
include qolr_lono/ldnsumbw1.txt
plot
include qolr_lono/ldnsumbw2.txt
plot
include qolr_lono/ldnsumbw3.txt
plot
include qolr_lono/ldnsumbw4.txt
plot
;----------------
pltcol=5 : select= : ydata=$5
include qolr_amod/ldnsumbw1.txt
plot
include qolr_amod/ldnsumbw2.txt
plot
include qolr_amod/ldnsumbw3.txt
plot
include qolr_amod/ldnsumbw4.txt
plot
;---------------- key ---------------------
msgdat=n : mhal=0 : mval=0 : msgsiz=1 : msglwt=1.1 : msgang=0 : mhkey=2.5
pltcol=0 : mvkey=0.3 : pltcol=1 : pltlwt=0.7 : mvsp=1.5
1.9 0.7 "
|11,pltcol=1||_4||11|  Leibold 2007
|11,pltcol=2||_4||11|  present study
"
mhkey=3.5 : pltlwt=0.9
2.4 1.8 "
|_0,pltcol=4| lono noise
|_0,pltcol=1| five tone
|_0,pltcol=2| flat noise
|_0,pltcol=5| AM tone
"
