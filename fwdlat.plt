; fwdlat.plt
;
clip=1 : head=0 : sizfac=5
XCYCLE  =     2.000
XMIN    =   250.000
XMAX    =  8000.000
XINT    =     0.   
YCYCLE  =     1.000
YMIN    =     1.000
YMAX    =    20.000
YINT    =     0.   
GRID    =     0.   
XLENGTH =     5.000
YLENGTH =     5.000
PLTYPE  =     1.000
TICDIR  =     1.000
SYMSIZ  =     1.000
MSGSIZ  =     -.13
LABSIZ  =     -.18
ANNSIZ  =     -.15
ANNLWT  =     1.000
LABLWT  =     1.000
PLTLWT  =     1.000
AXLWT   =     1.000
MSGLWT  =     1.000
YHOR    =     1.000
XANSKP  =     0.   
YANSKP  =     0.   
XPERCENT=    80.000
YPERCENT=    90.000
XLABEL  =frequency (Hz)
YLABEL  =forward latency (msec)
TLABEL  =auditory brainstem response
YANNOT  = 1 2 4 8 10 20
XFMT    =I4  
YFMT    =I4  
NEWFRAME
XLLC    =     1.750
YLLC    =     1.750

PLTYPE =     2.000
LINTYP=     0.   
pltcol=1 : lintyp=1 : pltlwt=1	
;------------------------------------
;include 1.fit
  ; y = a + b * c**(-i/100.) * f**(-d)
  ; a =  0.  , b =   13.0000  , c =   5.00000  , d =  0.410000
  ; f1 =  0.200000  , f2 =   8.00000   kHz
   200.000    18.227
  8000.000     4.017
  plot
   200.000    13.211
  8000.000     2.911
  plot
   200.000     9.575
  8000.000     2.110
  plot
   200.000     6.940
  8000.000     1.529
  plot
;------------------------------------
pltcol=1 : pltyp=3 : lintyp=0 : xdata=$1*1000 : zdata=0 : pltlwt=1
symb=2 : ydata=$6
include tbabr.txt
plot
pltcol=2
symb=1 : ydata=$7
include tbabr.txt
plot
pltcol=4
symb=4 : ydata=$8
include tbabr.txt
plot
pltcol=3
symb=5 : ydata=$9
include tbabr.txt
plot

msgsiz=1.1 : msgdat=0
3.0 4.6 "
|2,pltcol=1|   20
|1,pltcol=2|   40
|4,pltcol=4|   60
|5,pltcol=3|   80 dB SPL
"
0.3 0.6 "
|_1,pltcol=0| Neely et al. (1988)
"