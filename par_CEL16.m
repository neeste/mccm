function pa=par_CEL16
pa.parlab='CEL16';
pa.chsz=[1 1];                     % chamber size
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
pa.isv=[1136 1005 840 655 466 273 80]; % BM locations to save
pa.xtap=0.0; pa.xtex=6;pa.mmeq=1;
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
pa.ihceq=0; % 4 = Neely synapse model
pa.hbnl=0; pa.hbsc=0.04; pa.hbmx=5e-9; pa.ihcv=1;
pa.ihctc=0.2e-3; pa.ihc_s=1e5; pa.nrgn=2.5e6;
pa.ihc_A=14; pa.ihc_B=1; pa.ihc_C=17.836; pa.ihc_D=1.2637;
%---------- partition impedance
pa.k1o=2.358e+08; pa.k2o=3.003e+08; pa.k3o=3.154e+08; pa.k4o=4.112e+08;
pa.r1o=8.739e+02; pa.r2o=1.928e+03; pa.r3o=1; pa.r4o=0;
pa.m1o=7.257e-03; pa.m2o=3.286e-02;
pa.k1e=-2.733; pa.k2e=-3.350; pa.k3e=-2.942; pa.k4e=-2.833;
pa.r1e=-1.305; pa.r2e=-1.259; pa.r3e= 0.1; pa.r4e= 0.000;
pa.m1e=-0.061; pa.m2e=-0.082;
pa.xtap=0.2339; pa.hbt=0;
pa.aco=0.01; pa.ace=-0.4; pa.acq=0;
pa.bwo=0.05; pa.bwe=0; pa.bwq=0;
pa.gpo=1.00; pa.gpe=0; pa.gpq=0;
%--------------- partition defaults
pa.k1q=0; pa.r1q=0; pa.m1q=0;
pa.k2q=0; pa.r2q=0; pa.m2q=0; 
pa.k3q=0; pa.r3q=0;
pa.k4q=0; pa.r4q=0;
%---------- partition optmized
pa.ihc_A=10.784; pa.ihc_C=17.524; pa.ihc_D=1.9057;
pa.k1o=2.394e+08; pa.k2o=3.036e+08; pa.k3o=3.151e+08; pa.k4o=4.045e+08; 
pa.r1o=8.714e+02; pa.r2o=1.979e+03; pa.r3o=1.049e+00;
pa.m1o=7.417e-03; pa.m2o=3.417e-02;
pa.k1e=-2.6655; pa.k2e=-3.4762; pa.k3e=-2.9092; pa.k4e=-2.7948; 
pa.r1e=-1.3937; pa.r2e=-1.2466; pa.r3e=0.1033; 
pa.m1e=-0.0628*3; pa.m2e=-0.0828;
end % return
