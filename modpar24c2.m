function pa=modpar24c2
pa.chsz=[1 1];                 % two-chamber sizes
% CEL23 parameters --------------------------------------------
pa.parlab='cel24c2';
pa.gam = 1;                    % NDR multiplier
pa.n = 1401;                   % number of points along BM
pa.xl = 3.5;                   % scala length
pa.yw = 0.1;                   % scala width
pa.zh = 0.1;                   % scala height
pa.rho = 1;                    % fluid density
pa.bgn = 1.1;                  % B gain
pa.qgn=6.7e3;                  % qst gain
pa.stgain=0.0167;              % stimulus gain (target 40 dB SPL tone)
pa.isv=[1185 1046 838 657 489 307 121]; % BM locations to save
pa.xtap=0.0; pa.xtex=6;pa.mmeq=1;
pa.khe=0.0001; pa.rhe=0.0001; pa.mhe=0.01; % err=23.21 23.21
% OHC nolinearity parameters --------
pa.hbnl=0; pa.hbsc=0.05; pa.hbmx=1e-7;
% middle-ear parameters -------------
pa.nmev=4; % number [4=complete]
pa.kme=0.1; pa.rme=400; pa.mme=0.1; % simple ME
pa.acp=0.5; pa.adi=0.2; pa.aed=0.5; pa.ama=0.5; pa.ast=0.03;
pa.cep=2700;  pa.rep=29.16; pa.mevgn=1.262; pa.stim=3;
pa.mcp=0.0002; pa.rcp=0.1; pa.kcp=2200; pa.rfz=0.01;
pa.mdi=0.005;  pa.rdi=100;  pa.kdi=7.7e+06;
pa.mrw=5e-3; pa.rrw=20; pa.krw=1.5e5; pa.arw=0.0625;
pa.mco=3; pa.rco=1e5;  pa.gme=0.5; pa.vegn=10;
pa.kma=3e+06;  pa.rma=8000;  pa.mma=0.85;
pa.ked=7.04632e+06;  pa.red=200;  pa.med=0.02;
pa.kst=4.66022e+06;  pa.rst=8140.1;  pa.mst=0.89833;
pa.kim=3e+08;  pa.rim=28508.7; % nch=0 nme=4
% ---- parfit -------------------------------------------------
pa.k1o=6.681e+07;
pa.r1o=43.3864;
pa.m1o=0.0010604;
pa.k2o=3.3051e+07;
pa.r2o=9.23662;
pa.m2o=0.00290273;
pa.k3o=7.39901e+07;
pa.r3o=0.0917773;
pa.k4o=1.68115e+07;
pa.r4o=0;
pa.gpo=0.992953;
pa.aco=0.0219974;
pa.bwo=0.158727;
pa.k1e=-1.2195;
pa.r1e=0.0087;
pa.m1e=0.9975;
pa.k2e=-2.5197;
pa.r2e=2.4802;
pa.m2e=-0.0094;
pa.k3e=-5.9338;
pa.r3e=1.7912;
pa.k4e=0.0040;
pa.r4e=0.0000;
pa.gpe=-0.0029;
pa.ace=-0.4875;
pa.bwe=0.1836;
pa.k1q=-0.317379;
pa.r1q=-0.376951;
pa.m1q=-0.373343;
pa.k2q=-0.155603;
pa.r2q=-0.935460;
pa.m2q=0.278592;
pa.k3q=1.158179;
pa.r3q=-0.056186;
pa.k4q=-0.713082;
pa.r4q=0.000000;
pa.gpq=-0.000921;
pa.acq=0.000074;
pa.bwq=0.000089;
pa.hbt=0.0000;
end
