function pa=par_CEL22
pa.chsz=1;                     % chamber size
pa.parlab='CEL22';
pa.gam = 1;                    % NDR multiplier
pa.m = 1;                      % number of points across fluid
pa.n = 1401;                   % number of points along BM
pa.xl = 3.5;                   % scala length
pa.yw = 0.1;                   % scala width
pa.zh = 0.1;                   % scala height
pa.rho = 1;                    % fluid density
pa.bgn = 1;                    % B gain
pa.qgn=1;                      % qst gain
pa.stgain=0.0167;              % stimulus gain (target 40 dB SPL tone)
pa.isv = [1200 1019 822 515 80 40 20]; % BM locations to save
pa.xtap=0.0; pa.xtex=6;pa.mmeq=1;
pa.khe=0.0001; pa.rhe=0.0001; pa.mhe=0.01; % err=23.21 23.21
% middle-ear parameters -------------
pa.nmev=4; % number [4=complete]
pa.kme=0.1; pa.rme=400; pa.mme=0.1; % simple ME
pa.acp=0.5; pa.adi=0.2; pa.aed=0.5; pa.ama=0.5; pa.ast=0.03*100;
pa.cep=2700;  pa.rep=29.16; pa.mevgn=1.262; pa.stim=3;
pa.mcp=0.0002; pa.rcp=0.1; pa.kcp=2200; pa.rfz=0.01;
pa.mdi=0.005;  pa.rdi=100;  pa.kdi=7.7e+06;
pa.mrw=5e-3; pa.rrw=20; pa.krw=1.5e5; pa.arw=0.0625;
pa.mco=3; pa.rco=1e5;  pa.gme=0.5; pa.vegn=100;
pa.kma=3e+06;  pa.rma=8000;  pa.mma=0.85;
pa.ked=7.04632e+06;  pa.red=200;  pa.med=0.02;
pa.kst=4.66022e+06;  pa.rst=8140.1;  pa.mst=0.89833;
pa.kim=3e+08;  pa.rim=28508.7; % nch=0 nme=4
% ---- parfit ----
pa.k1o=1.054e+08;
pa.r1o=113.1;
pa.m1o=0.002;
pa.k2o=9.443e+06;
pa.r2o=15.99;
pa.m2o=0.002;
pa.k3o=1.767e+07;
pa.r3o=0.831;
pa.k4o=2.482e+07;
pa.r4o=0;
pa.aco=0.01;
pa.k1e=-0.9662;
pa.r1e=0.8793;
pa.m1e=1.7012;
pa.k2e=-2.0918;
pa.r2e=1.0379;
pa.m2e=-0.2339;
pa.k3e=-3.9166;
pa.r3e=2.5927;
pa.k4e=0.2348;
pa.r4e=0.0000;
pa.ace=0.0000;
pa.k1q=-0.290055;
pa.r1q=-0.320422;
pa.m1q=-0.583975;
pa.k2q=-0.503048;
pa.r2q=-0.540285;
pa.m2q=0.646131;
pa.k3q=0.831386;
pa.r3q=-0.359765;
pa.k4q=-0.598210;
pa.r4q=0.000000;
pa.acq=0.000000;
%---------------
pa.bwo=0.05; pa.bwe=0; pa.bwq=0;
pa.gpo=1.00; pa.gpe=0; pa.gpq=0;
end % return
