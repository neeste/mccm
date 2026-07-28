; amtmsk.plt
;
head=0 : ticsiz=0.7 : ticdir=in : yhor=y
axlwt=0.7 : annlwt=1.1 : lablwt=1.1 : pltlwt=0.9 : labsiz=1.1 : sizfac=6
pltyp=line : clip=y
xdata=$1 : zdata=0
;
newframe
xllc=1.5 : yllc=2.0
xlen=3.5 : ylen=3.5
xlabel=masker level (dB SPL)
ylabel=amount of masking (dB)
tlabel=forward masking
xanskp=0 : yanskp=0
pltyp=lines : pltlwt=1
xmin=40 : xmax=80 : xint=4.2 : xper=80
ymin=0  : ymax=40 : yint=8   : yper=100
;
pltyp=lines : lintyp=6 : pltlwt=0.7; pltcol=0
ydata=$2
include amtmsk0.txt
plot 
ydata=$3
include amtmsk0.txt
plot
ydat=$4
include amtmsk0.txt
plot
;
pltyp=both : lintyp=0 : pltlwt=1
pltcol=1 : ydata=$2
include amtmsk.txt
plot 
pltcol=2 : ydata=$3
include amtmsk.txt
plot
pltcol=4 : ydat=$4
include amtmsk.txt
plot
;
msgdat=0 : pltlwt=1
0.5 3.4 "
probe delay
|_0,pltcol=1|  10
|_0,pltcol=2|  20
|_0,pltcol=4|  40
"
1.2 0.3 "|_6,pltcol=0,pltlwt=0.7| Jesteadt"

