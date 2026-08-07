% AMPWHY -- WHY is rlsplit's amplifier dead? A mechanism, not another observation.
%
% HYPOTHESIS. The active force is proportional to HAIR BUNDLE displacement, and
% the hair bundle is driven by the TM/RL shear d2. rlsplit sets Df(2,:) = 0, i.e.
% it takes the fluid drive AWAY from d2 and gives it to d3. If d2 then barely
% moves, the active force has almost nothing to act on and the amplifier is dead
% for a STRUCTURAL reason, which would explain why 300 evaluations with the gain
% shortfall as the leading term moved it only 8.21 -> 11.42 dB.
%
% PREDICTION: d2mx collapses under rlsplit relative to d1mx, and ohcBM (net OHC
% work into the BM, >0 = amplifying) collapses with it.
% REFUTATION: d2mx holds up, in which case the dead amplifier is something else
% and the BM-RL force-pair idea below is not supported.
PROJ = '/Users/neely/Library/CloudStorage/OneDrive-FatherFlanagan''sBoysHome/STN/mod/mccm';
addpath(PROJ); cd(PROJ);
V = {'2-DOF (arm B final)', '/Users/neely/mccm_runs/b2fit_seg2.mat', struct('d3int',0,'rlsplit',0); ...
     'd3int internal',      '/Users/neely/mccm_runs/d3stage_B_seg2.mat', struct('d3int',1,'rlsplit',0); ...
     'rlsplit fitted',      '/Users/neely/mccm_runs/rlfit_seg2.mat', struct('d3int',1,'rlsplit',1)};
fprintf('\n  %-20s %10s %10s %10s %12s %10s\n','config','d1mx','d2mx','d2/d1','ohcBM','amp dB');
for i=1:size(V,1)
    L=load(V{i,2}); pa=L.R.pa;
    if (isfield(pa,'hbmode')), pa=rmfield(pa,'hbmode'); end
    f=fieldnames(V{i,3}); for j=1:numel(f), pa.(f{j})=V{i,3}.(f{j}); end
    try
        S=tdm26('wnr1',struct('fr',2,'lv',60,'pa',pa),0,0); d=S.dgn;
        s=score26(pa,'fast',false);
        fprintf('  %-20s %10.3e %10.3e %10.4f %12.3e %+10.2f\n', V{i,1}, ...
                d.d1mx, d.d2mx, d.d2mx/max(d.d1mx,eps), d.ohcBM, s.amp_gain);
    catch e, fprintf('  %-20s FAILED %s\n', V{i,1}, e.message(1:min(50,end))); end
end
fprintf('\n  d2/d1 is the shear per unit BM motion, i.e. how much drive the hair\n');
fprintf('  bundle gets. ohcBM > 0 means the OHC delivers energy TO the BM.\n');
disp('AMPWHY_DONE');
