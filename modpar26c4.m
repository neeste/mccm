function pa=modpar26c4
pa.parlab='cel26c4';
% ---- NESTED CONSTRUCTION, default as of 2026-07-26 -----------------------
% CL is CARVED FROM SV in BOTH senses, which is one decision and not two:
%   AREA      SV gives up CL's share, so the total is CONSERVED at 2.00. The
%             old [0.95 0.05 1.0 0.05] summed to 2.05, i.e. the 4th chamber was
%             APPENDED and grew the cross-section. The cochlear cross-section is
%             fixed, so resolving CL must subdivide it.
%   TOPOLOGY  CL takes the position next to the BM, so d1 spans ST-CL, not
%             ST-SV. Stack: ST | BM(d1) | CL | RL(d3) | SS --d2-- SV.
pa.chsz=[0.95 0.05 0.95 0.05]; % ST SS SV CL  (CL carved from SV, total 2.00)
pa.nested = 1;                 % d1 spans ST-CL (see above)
pa.clvtgt = 2;                 % vent CL -> SS. NOT SV and emphatically NOT ST:
                               % a vent stamped onto entries already carrying a
                               % partition stiffens that partition's own
                               % pressure difference instead of opening a path.
                               % ST (entries 4,13 = d1) DESTROYS the model
                               % (maxRe 23263). SV works but is 5x worse.
pa.clvent = 3.0;               % vent strength; sets channel inertance m5/clvent
pa.clvoct = 4.0;               % CL resonance, octaves below local CF. Resolved
                               % per place in tdm26/fdm26 from k1/m1/m5, so it
                               % tracks CF and survives a refit of those.
% MEASURED (score26 fast, unfitted -- these carry m=3b's parameter set):
%                             maperr    amp d1    maxRe  range  fold   chi
%   appended, the old default 1020.3   +80.62  +4524.8   0.91  0.91  12.1
%   nested + SS vent           522.7   +55.08     +3.9   4.52  1.43   4.2
%   + clvoct 4 (THIS DEFAULT)  525.2   +56.59     +2.3   4.44  1.28   6.4
%   m=3b scaffolding (FITTED)  499.3   +81.15    +19.3   5.74  0.84   9.2
% All re-measured 2026-07-27 with the corrected BF detector and cross-checked by
% two independent runs. The appended row's chi 12.1 is meaningless: that model
% has a collapsed map (range 0.91) and maxRe +4525.
%
% BF DETECTOR, for reading the range/fold columns. An apical EDGE GUARD was
% added: the old detector took a global argmax of an f-weighted curve, so at
% apical places, where the response above CF is negligible, it returned the band
% edge -- two places both reporting 24.98 kHz at -39 dB. It now returns NaN
% there, which discards ~25% of places (all apical), so RANGE UNDERSTATES the
% true span. A continuity-constrained tracker was also tried and REVERTED: it
% followed the sub-dominant mode past a crossing, to rank 11 and -16 dB. See the
% long note in score26/local_tip before re-attempting it.
% FOLD MEASURES MODE DEGENERACY, NOT MAP QUALITY. m=3b's 0.84 oct is the BM mode
% crossing the shear mode at x/L 0.10 with the runner-up 0.01 dB down. Doubling
% SS to [0.95 0.10 1.00] uncrosses them (fold 0.09, a clean map) and improves
% maperr to 406.8 -- see bisect_area.m; that is a 3-CHAMBER result and is NOT
% applied here.
% clvoct trades tip-tail contrast against the CF map and 4.0 is the free point:
% chi 4.2 -> 6.4, a 52% gain, for a 0.5% maperr cost (522.7 -> 525.2). Both
% halves re-verified with the corrected detector.
% 2.0 oct gives markedly more contrast but costs 55% maperr (522.7 -> 811.0);
% 1.0 oct puts the map back at appended levels. SN's 0.5 oct prescription leaves
% no amplifier at all (amp +0.05, maxRe +145.2) -- the mechanism is real but the
% model's usable placement is much further from CF than the neural data suggest,
% and that gap is unexplained.
% To recover the previous behaviour: pa.nested=0, and remove clvent/clvoct.
% Saved fits loaded from .mat carry none of these fields and so still run the
% appended construction, meaning nothing already banked changes.
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
%---------- DOF-3 (OC height / cortilymph pump) impedance  [PROVISIONAL]
% d3 = change in BM-RL separation (organ-of-Corti height), driven by OHC somatic
% motility; couples CL<->SS.  Seeded from the DOF-2 profile as a neutral,
% physically-scaled starting point -- these are the first parameters to fit.
pa.k5o=3.05084e+08; pa.k5e=-3.5063; pa.k5q= 0.000002;
pa.r5o=    1984.19; pa.r5e=-1.2479; pa.r5q= 0.000002;
pa.m5o=  0.0360276; pa.m5e=-0.0812; pa.m5q= 0.000002;
end % return
