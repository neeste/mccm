% Do the multi-chamber configs share the m=1 champion's (better-tuned) impedance
% profile, or do they start from separate defaults? Cheap .mat load, no compute.
L=load('refit_c1broad.mat'); p1=L.R.pa;
p3=modpar26(3); p4=modpar26(4);
L3=load('refit_c3_map.mat'); f3=L3.R.pa;
nm={'k1o','k1e','k1q','m1o','m1e','r1o','r1e','k3o','k3e','k4o','k4e','gpo','gpe','aco','ace','gam'};
fprintf('\n  param |    m=1 champ |   m=3 fitted |  m=4 default | c1 vs c4 ratio\n');
fprintf('%s\n',repmat('-',1,72));
for i=1:numel(nm)
    n=nm{i};
    v1=NaN; v3=NaN; v4=NaN;
    if isfield(p1,n), v1=p1.(n); end
    if isfield(f3,n), v3=f3.(n); end
    if isfield(p4,n), v4=p4.(n); end
    r=NaN; if (isfinite(v1)&&isfinite(v4)&&v4~=0), r=v1/v4; end
    fprintf('  %-5s | %12.5g | %12.5g | %12.5g | %8.3f\n', n, v1, v3, v4, r);
end
fprintf('\n(ratio far from 1.0 = the m=1 champion sits somewhere the m=4 default does not)\n');
disp('SHARED_PAR_DONE');
