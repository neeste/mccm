function pa=modpar26(nch)
% chamber sizes
if (nch<3)
    pa=par_CEL16;
elseif (nch==3)
    pa=modpar26c3;
    % --- m=3 CONSOLIDATED with the internal third DOF, on by default ---------
    % d3 here is INTERNAL: no fluid compartment, and NO BACK-ACTION (s1 and s2
    % read only d1,d2,v1,v2; d3 enters no fluid equation; a(i1)/a(i2) never see
    % it). Verified bit-identical to the 2-DOF m=3 on every scored quantity --
    % maperr, maxRe, maxRe_osc, map, contrast and amplifier all matched to zero
    % -- so this is a strict superset, not a model change. Cost is 27%.
    % What it buys: d3 is the absolute RETICULAR LAMINA position, available on
    % every 3-chamber run. Shadow prediction RL/BM ~ 1.52, which is the
    % coordinate the mouse OCT work concerns.
    % Set pa.d3int=0 after the call for the legacy 2-DOF behaviour. Saved fits
    % loaded from .mat carry no d3int field and therefore run 2-DOF, producing
    % identical results either way, so nothing already banked changes.
    % k5/r5/m5 seeded from modpar26c4 (provisional there); they are the first
    % things to fit once d3 acquires a compartment at m=4.
    pa.d3int = 1;
    pa.k5o = 3.05084e+08; pa.k5e = -3.5063; pa.k5q = 0.000002;
    pa.r5o =     1984.19; pa.r5e = -1.2479; pa.r5q = 0.000002;
    pa.m5o =   0.0360276; pa.m5e = -0.0812; pa.m5q = 0.000002;
else
    pa=modpar26c4;   % 4th chamber = CL (cortilymph space), 3rd partition DOF
end
pa.m=nch;             % number of fluid chambers
pa.xp=1;              % length of excitation pattern
% TDM parameters
pa.nstp = 20480; pa.ncyc = 1; pa.ntsw = 10; pa.ntsf = 2048; pa.dt = 2e-6; pa.nimp = 1;
pa.ME_VEP_SCALE = 10; pa.ME_VST_SCALE_PHANTOM = 2; pa.IHC_SUBSTEPS = 10;
pa.aflom_fac = 1;     % fluid-coupling divisor factor (ecochg pins to 2 for tdm25 parity)
pa.gampro = ones(pa.n,1); 
pa.synpro = ones(pa.n,1);
pa.met=200; pa.hco=1e-4; pa.hcs=1e-4;
end

%--------------------------------------------------------------

function pa=par_CEL16
pa.parlab='CEL16';
pa.chsz=[1 1];                 % chamber size
pa.gam = 1;                    % NDR multiplier
pa.m = 1;                      % number of points across fluid
pa.n = 1401;                   % number of points along BM
pa.xl = 3.5;                   % scala length
pa.yw = 0.1;                   % scala width
pa.zh = 0.1;                   % scala height
pa.rho = 1;                    % fluid density
pa.stgain=0.0127;              % stimulus gain (target 40 dB SPL tone)
pa.isv=[1136 1005 840 655 466 273 80]; % BM locations to save
pa.hbt=8.5; pa.xtap=0.2339; pa.xtex=6; pa.mmeq=1;
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
pa.ldew=2; pa.ldne=140; pa.ldsc=[1.360 7.057 3.518 18.93 63.79];
%---------- partition impedance
pa.k1o=  2.394e+08; pa.k1e=-2.6655; pa.k1q= 0.000000; % CEL16
pa.r1o=      871.4; pa.r1e=-1.3937; pa.r1q= 0.000000;
pa.m1o=   0.007417; pa.m1e=-0.1884; pa.m1q= 0.000000;
pa.k2o=  3.036e+08; pa.k2e=-3.4762; pa.k2q= 0.000000;
pa.r2o=       1979; pa.r2e=-1.2466; pa.r2q= 0.000000;
pa.m2o=    0.03417; pa.m2e=-0.0828; pa.m2q= 0.000000;
pa.k3o=  3.151e+08; pa.k3e=-2.9092; pa.k3q= 0.000000;
pa.r3o=      1.049; pa.r3e= 0.1033; pa.r3q= 0.000000;
pa.k4o=  4.045e+08; pa.k4e=-2.7948; pa.k4q= 0.000000;
pa.r4o=          0; pa.r4e= 0.0000; pa.r4q= 0.000000;
pa.gpo=          1; pa.gpe= 0.0000; pa.gpq= 0.000000;
pa.aco=       0.01; pa.ace=-0.4000; pa.acq= 0.000000;
pa.bwo=       0.05; pa.bwe= 0.0000; pa.bwq= 0.000000;
end % return
