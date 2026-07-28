% ldngrw - loudness growth protocol set
st=1; % stimulus type [0=multi-tone 1=flat-noise 2=lono-noise 3=AM_tone]
tdm24(0,1,0,0)       % Ped calibration
tdm24('ldngrw')      % single tone only
tdm24('ldngrw',1,st) % loudness summation
tdm24('ldngrw',2,st) % ensemble width
tdm24('ldngrw',3,st) % number of ensemble
tdm24('ldngrw',4,st) % loudness summation @ 60 dB SPL
