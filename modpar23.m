function pa=modpar23(nc)
% chamber sizes
if (nc==2)
    pa=modpar23c2;
elseif (nc==3)
    pa=modpar23c3;
end
end

function pa=modpar23c2
% chamber sizes
pa.chsz=[1 1];
% default parameters
pa.gam= 1.000;
pa.n=2501;
pa.xl = 2.500;
pa.yw = 0.100;
pa.zh = 0.100;
pa.rho= 1.000;
pa.bwo= 0.050;
pa.bwe= 0.000;
pa.hbt= 20; pa.xp=0.97;
pa.xtap=0.0; pa.xtex=6;
pa.ast= 0.01;   % area of stapes
pa.ahe= 0;      % area of helicotrema
pa.khe=1.525; pa.rhe=1.096e-09; pa.mhe=55.89928;
pa.kme=5.119e+06; pa.rme=528.54; pa.mme=0.061404; % err=43.81 43.13
% ---- parfit ----
pa.k1o=3.724e+05;
pa.r1o=0.02581;
pa.m1o=2.025e-06;
pa.k2o=3.391e+06;
pa.r2o=8.975;
pa.m2o=3.987e-05;
pa.k3o=5.076e+05;
pa.r3o=5.396e-05;
pa.k4o=5.671e+05;
pa.r4o=2e-06;
pa.aco=0.01;
pa.k1e=-2.4274;
pa.r1e=3.1694;
pa.m1e=0.0013;
pa.k2e=-6.5431;
pa.r2e=-5.2047;
pa.m2e=-2.2354;
pa.k3e=-4.0360;
pa.r3e=2.0002;
pa.k4e=-3.8753;
pa.r4e=0.0000;
pa.ace=6.7599;
pa.k1q=0.191319;
pa.r1q=-1.273363;
pa.m1q=-0.000179;
pa.k2q=0.000019;
pa.r2q=1.530981;
pa.m2q=0.657814;
pa.k3q=-0.385552;
pa.r3q=0.591947;
pa.k4q=0.831565;
pa.r4q=0.000000;
pa.acq=-2.101424;
end

function pa=modpar24c3
pa.parlab='cel24c3';
pa.chsz=[0.9 0.1 1];           % three-chamber sizes
pa.gam = 1;                    % NDR multiplier
pa.m = 1;                      % number of points across fluid
pa.n = 1401;                   % number of points along BM
pa.xl = 3.5;                   % scala length
pa.yw = 0.1;                   % scala width
pa.zh = 0.1;                   % scala height
pa.rho = 1;                    % fluid density
pa.stgain=0.0127;              % stimulus gain (target 40 dB SPL tone)
pa.isv=[1136 1005 840 655 466 273 80]; % BM locations to save
pa.xtap=0.0; pa.xtex=6;pa.mmeq=1;
pa.bwo = 0.05; pa.bwe = 0; pa.bwq = 0; % BM width
pa.gpo = 1.00; pa.gpe = 0; pa.gpq = 0; % partition gain (HB re BM)
pa.hbt=0; pa.xp=1;
pa.khe=0.0001; pa.rhe=0.0001; pa.mhe=0.01; % err=23.21 23.21
% middle-ear parameters -------------
pa.nmev=4; pa.kme=0.1; pa.rme=400; pa.mme=0.1;
pa.acp=0.5; pa.adi=0.2; pa.aed=0.5; pa.ama=0.5; pa.ast=0.03;
pa.cep=2700;  pa.rep=29.16; pa.mevgn=1.262; pa.stim=3;
pa.mcp=0.0002; pa.rcp=0.1; pa.kcp=2200; pa.rfz=0.01;
pa.mdi=0.005;  pa.rdi=100;  pa.kdi=7.7e+06;
pa.mrw=5e-3; pa.rrw=20; pa.krw=1.5e5; pa.arw=0.0625;
pa.mco=3; pa.rco=1e5;  pa.gme=0.5;
pa.kma=3e+06;  pa.rma=8000;  pa.mma=0.85;
pa.ked=7.04632e+06;  pa.red=200;  pa.med=0.02;
pa.kst=4.66022e+06;  pa.rst=8140.1;  pa.mst=0.89833;
pa.kim=3e+08;  pa.rim=28508.7; % nch=0 nme=4
% hair-cell & neural parameters --------
pa.ihceq=4; % 4 = Neely synapse model
pa.hbnl=0; pa.hbsc=0.04; pa.hbmx=6e-9; pa.ihcv=1;
pa.ihctc=0.2e-3; pa.ihcsf=1e5; pa.nrgn=2.5e6; pa.ew=0.1;
pa.ihcex=14; pa.ihcrr=17.836; pa.ihcdr=1.2637;
% ---- parfit ----
pa.k1o=3.765e+05;
pa.r1o=0.0273;
pa.m1o=2.026e-06;
pa.k2o=3.384e+06;
pa.r2o=7.846;
pa.m2o=3.873e-05;
pa.k3o=4.885e+05;
pa.r3o=5.649e-05;
pa.k4o=4.708e+05*0.5;
pa.r4o=2e-06;
pa.aco=0.01;
pa.k1e=-2.4410;
pa.r1e=3.2939;
pa.m1e=0.0014;
pa.k2e=-6.4479;
pa.r2e=-5.2853;
pa.m2e=-2.2830;
pa.k3e=-4.0139;
pa.r3e=1.9891;
pa.k4e=-4.0853;
pa.r4e=0.0000;
pa.ace=6.5572;
pa.k1q=0.172199;
pa.r1q=-1.199597;
pa.m1q=-0.000223;
pa.k2q=0.000025;
pa.r2q=1.535546;
pa.m2q=0.714623;
pa.k3q=-0.343642;
pa.r3q=0.587935;
pa.k4q=0.945069;
pa.r4q=0.000000;
pa.acq=-2.177005;
end
