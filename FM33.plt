; fm2.plt - function fit to loudness vs level data from FM33, Table III

yhor=y : head=0 : clip=y : ticdir=in : sizfac=7

pltyp=line : grid=y
xmin=-10 : xmax=130 : xlen=7 : xint=14 : xanskp=1.1
ymin=0.01 : ymax=1000000 : ylen=6 : ycyc=1 : yper=100
xllc=2 : yllc=1
tlabel=Fletcher and Munson (1933, Table III)
xdat=$1 : ydat=$2
xlab=intensity (dB SPL)
ylab=loudness (millisones)
include FM33.txt
plot


