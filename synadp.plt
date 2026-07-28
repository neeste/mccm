; fwdmsk.plt - forward masking: masker and probe
;
head=0 : ticsiz=0.7 : ticdir=in : yhor=y
axlwt=0.7 : annlwt=1.1 : lablwt=1.1 : pltlwt=0.9 : labsiz=1.1 : sizfac=6
pltyp=line : clip=y
xdata=$1 : zdata=0
lintyp=0 : ydata=$2/1e10
;
newframe
xllc=1.5 : yllc=2.0
xlen=3.5 : ylen=3.5
xanskp=0
xmin=0 : xmax=120 : xint=16 : xper=100*8/9 : xanskp=3
ymin=0 : ymax=3 : yint=3.2 : yper=75
xlabel=time (msec)
ylabel=WNR (spks/sec)
tlabel=
pltcol=1
include fwdmsk32.txt
plot
pltcol=2
include fwdmsk22.txt
plot
pltcol=4
include fwdmsk12.txt
plot
;---------------- key ---------------------
msgdat=n : mhal=0 : mval=0 : msgsiz=1 : msglwt=1.1 : msgang=0 : mhkey=3.3
1.8 3.2 "
|_0,pltcol=1|  80
|_0,pltcol=2|  60
|_0,pltcol=4|  40
"
0.1 3.6 "x\ 10{10}"
