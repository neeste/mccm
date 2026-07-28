; test_band.plt - quarter-octave stimuli
;
head=0 : ticsiz=0.7 : ticdir=in : yhor=y
axlwt=0.7 : annlwt=1.1 : lablwt=1.1 : pltlwt=0.9 : labsiz=1.1 : sizfac=6
pltyp=line : clip=y
xdata=$1 : ydata=$2 : zdata=0 : lintyp=0
xllc=1 : yllc=0
xlen=6 : ylen=1.5
xanskp=0
xmin=0 : xmax=100 : xint=10 : xper=100
ymin=-3 : ymax=3 : yint=6 : yper=60 : yshif=-10
yann=-3 0 3
xlabel=time (msec)
tlabel=
tmin=0 : tmax=100 : tint=1 : tanskp=-1
msgdat=y;--------------------------
newframe
pltcol=5 : yllc=(yllc+ylen)
include test_band_4.txt ; AM tone
plot10 4 "AM tone"
;--------------------------
xanskp=-1 : xlabel=
;--------------------------
newframe
pltcol=2 : yllc=(yllc+ylen)
include test_band_2.txt
plot
10 4 "flat noise"
;--------------------------
newframe
pltcol=1 : yllc=(yllc+ylen)
include test_band_1.txt
plot
10 4 "five-tone complex"
;--------------------------
newframe
pltcol=4 : yllc=(yllc+ylen)
include test_band_3.txt ; low-noise noise
plot
10 4 "low-noise noise"
;--------------------------
