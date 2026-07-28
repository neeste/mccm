% Proper 3-chamber joint fit: open the structural lever (chsz) + BM/active,
% warm-started near the chsz2 value where the SS-size sweep crossed the target
% slope (~0.37), so fminsearch refines the multi-objective balance there.
R.pa = modpar26(3); R.pa.chsz = [0.95 0.37 1.0]; save('warm_c3.mat','R');
fprintf('warm-start chsz=[0.95 0.37 1.0]\n');
% fit: k1o(1) m1o(3) aco(10) k1e(11) m1e(13) ace(20) + chsz(31 32 33)
refit_tdm(3, [1,3,10,11,13,20,31,32,33], 150, 'refit_c3_proper.mat', 'warm_c3.mat');
