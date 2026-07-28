; PLT file created by Greg2000(2.4.0)
; Date and time created:    5/3/2005 9:26:01 AM
;
;============  PAGE    1  =============
NEWPAGE

HEADER=0
SIZFAC  =6.000
MAGNIFY =1.000
ROTATE  =0.000
VXLLC   =0.000
VYLLC   =0.000
XCYCLE = 1
XMIN = 0.03
XMAX = 3
xdata=$1/1000
XINT = 10
xann=0.03 0.1 0.3 1 3
YCYCLE = 0.000
YMIN = 50
YMAX = 70
YINT = 4
GRID    =0.000
XLENGTH =5.000
YLENGTH =5.000
PLTYPE = 1.000
TICDIR  =1.000
SYMSIZ  =1.200
MSGSIZ  =1.000
LABSIZ  =1.200
ANNSIZ  =1.200
ANNLWT  =1.000
LABLWT  =1.000
PLTLWT  =1.000
ECHO    =0.000
AXLWT   =1.000
MSGLWT  =1.000
YHOR    =1.000
XANSKP  =0.000
YANSKP  =0.000
XPERCENT=80.000
YPERCENT=100
XLABEL = 
YLABEL = 
XFMT = I4
YFMT = I4
TICSIZ  =1.000

NEWFRAME
XLLC    =1.750
YLLC    =1.750
XANSKP  =0.000
YANSKP  =0.000
MXLLC   =0.250
MYLLC   =2.200
MSGSIZ  =1.400
MHALIGN = 2
MSGANGLE = 90.000
-0.9 2.5 "level (dB SPL) of 1-kHz tone 
judged equally loud to complex"
MSGANGLE = 0.000
MXLLC = 2.500 
MYLLC = -0.700
"bandwidth of complex (kHz)"
MHALIGN = 0
SYMBOL = 11.000
LINTYPE = 0.000
zdata=0
PLTYPE = 2
DATA
    30.045     55.664
    141.457     55.664
PLOT

PLTYPE = 1
DATA
    46.000     56.071     1.21
    92.000     55.257     1.94
PLOT


PLTYPE = 2
DATA
    141.457      55.664
    2642.409     67.693
PLOT

PLTYPE = 1
DATA
    231.000     57.829     1.53
    465.000     60.326     1.54
    956.000     63.529     1.79
    2119.000     66.850     .79
PLOT

PLTYPE = 2
LINTYPE = 1.000
DATA
    132.6        50
    132.6        65
PLOT

SYMBOL = 1.000
LINTYPE = 0.000

PLTYPE = 2
DATA
    30.045     55.664
    141.457     55.664
PLOT

PLTYPE = both
lintype = 5
symbol=0
zdata = 0
clip=y
pltlwt=1.2 : pltcol=1 : symbol=0
;include bwlmtc-nh.txt
;plot

pltcol=0
data
msgsiz=1.1
MXLLC = 0.5
MYLLC = 4.7
mvsp = 1.8
"
|11,pltlwt=1,pltcol=0| Leibold et al. data
"
