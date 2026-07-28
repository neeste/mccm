% Compute the basal-boundary impedance ratio for the 3-chamber tbabr config.
% At the base (k=1, x=0) the interior basal coupling L1_p is REPLACED by 2*alfx.
% Matched (non-reflecting) requires 2*alfx ~ the coupling it replaces (~L1_c).
pa=modpar26(3);
fprintf('nch(m)=%d  aflom_fac=%g  chsz(raw)=[%s]\n', pa.m, pa.aflom_fac, num2str(pa.chsz));
n=pa.n; dx=pa.xl/(n-1); rho=pa.rho; af=pa.aflom_fac;
bw1=pa.bwo; m11=pa.m1o; ac1=pa.aco;         % x=0 -> all exp(...) = 1
abmom1=bw1*dx/m11;
aflom1=ac1/(af*rho*dx);
chsz=pa.chsz*(2/sum(pa.chsz));
alfx =(pa.ast/pa.mst)/abmom1;
L1c  =aflom1*chsz(1)/abmom1;
fprintf('dx=%.4g  abmom(1)=%.4g  aflom(1)=%.4g  chsz(norm)=[%s]\n',dx,abmom1,aflom1,num2str(chsz));
fprintf('ast=%.4g mst=%.4g  ->  alfx=%.4g   2*alfx=%.4g\n',pa.ast,pa.mst,alfx,2*alfx);
fprintf('line coupling L1_c(1)=%.4g   (interior basal coupling L1_p ~ L1_c)\n',L1c);
fprintf('----\nBOUNDARY / LINE  2*alfx / L1_c = %.4g\n', 2*alfx/L1c);
fprintf('(=1 matched;  >>1 rigid/reflecting;  <<1 pressure-release/reflecting)\n');
disp('BOUNDARY_DIAG_DONE');
