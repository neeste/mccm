function pa=modpar26c3
pa.parlab='cel26c3';
pa.chsz=[0.95 0.05 1.0]; % err=93.36 93.36
pa.gam = 1;                    % NDR multiplier
pa.m = 1;                      % number of points across fluid
pa.n = 1401;                   % number of points along BM
pa.xl = 3.5;                   % scala length
pa.yw = 0.1;                   % scala width
pa.zh = 0.1;                   % scala height
pa.rho = 1;                    % fluid density
pa.stgain=0.0127;              % stimulus gain (target 40 dB SPL tone)
pa.isv=[1136 1005 840 655 466 273 80]; % BM locations to save
pa.hbt=33; pa.xtap=0.2339; pa.xtex=6; pa.mmeq=1;
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
pa.ldew=2; pa.ldne=140; pa.ldsc=[1.142 8.296 4.841  8.05 34.00];
%---------- partition impedance
pa.k1o=2.41442e+08; pa.k1e=-2.7377; pa.k1q= 0.000002;
pa.r1o=     1056.6; pa.r1e=-1.3138; pa.r1q= 0.000002;
pa.m1o= 0.00745053; pa.m1e=-0.2012; pa.m1q= 0.000002;
pa.k2o=3.05084e+08; pa.k2e=-3.5063; pa.k2q= 0.000002;
pa.r2o=    1984.19; pa.r2e=-1.2479; pa.r2q= 0.000002;
pa.m2o=  0.0360276; pa.m2e=-0.0812; pa.m2q= 0.000002;
pa.k3o=3.15082e+08; pa.k3e=-2.9231; pa.k3q= 0.000003;
pa.r3o=    1.10392; pa.r3e= 0.1065; pa.r3q= 0.000002;
pa.k4o=4.04629e+08; pa.k4e=-2.8017; pa.k4q= 0.000001;
pa.r4o=          0; pa.r4e= 0.0000; pa.r4q= 0.000000;
pa.gpo=     0.9956; pa.gpe= 0.0000; pa.gpq= 0.000003;
pa.aco= 0.00994045; pa.ace=-0.4324; pa.acq= 0.000002;
pa.bwo=  0.0517652; pa.bwe= 0.0001; pa.bwq= 0.000002;
end % return
