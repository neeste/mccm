function R = fit_dof3(maxfe)
%FIT_DOF3  Fit the third-DOF parameters k5/r5/m5 at m=3b.
%
% These have NEVER been fitted. They were structurally inert in BOTH solvers
% until 2026-07-29: exactly 0.000e+00 sensitivity, because d3 had no back-action
% (tdm26) and no third DOF at all (fdm26's m==3 branch). parfit26's parnames()
% lists 30 parameters and none of them is k5/r5/m5, so the fit machinery never
% exposed them either. m5o is still the placeholder I seeded from m2o.
%
% OBJECTIVE: fd maperr at m=3b, plus a STABILITY penalty. The stability term is
% not optional: score26 'fast' runs fdm26, which is frequency-domain and cannot
% see the time march diverge. That is exactly how m=4 stayed green for a session
% while tdm26 blew up. Checked on a coarse grid to keep the cost down.
if (nargin<1 || isempty(maxfe)), maxfe = 120; end
p0 = modpar26(3);                      % m=3b, third DOF live
nm = {'k5o','r5o','m5o','k5e','r5e','m5e'};
v0 = cellfun(@(f) p0.(f), nm);
% o-params are positive scales -> fit in log; e-params are exponents -> additive.
x0 = [0 0 0 v0(4) v0(5) v0(6)];
ISV = [1136 1005 840 655 466 273 80];

    function pa = setp(x)
        pa = p0;
        pa.k5o = v0(1)*exp(x(1)); pa.r5o = v0(2)*exp(x(2)); pa.m5o = v0(3)*exp(x(3));
        pa.k5e = x(4);            pa.r5e = x(5);            pa.m5e = x(6);
    end

    function J = obj(x)
        pa = setp(x);
        try
            S = score26(pa,'fast',false); J = S.maperr;
            if (~isfinite(J)), J = 1e6; return; end
        catch, J = 1e6; return; end
        % STABILITY GUARD, every 6th evaluation (tdm26 is ~13 s, fdm ~1 s).
        persistent cnt; if (isempty(cnt)), cnt = 0; end
        cnt = cnt + 1;
        if (mod(cnt,6)==0)
            q = pa; q.isv = ISV;
            try
                evalc('T = tdm26(0,q,0,0);'); d = T.d1(:);
                if (~all(isfinite(d)) || max(abs(d)) > 1e-1), J = J + 1e5; end
            catch, J = J + 1e5; end
        end
    end

J0 = obj(x0);
fprintf('\n  start J (maperr) = %.4f\n', J0);
fprintf('  fitting %s\n', strjoin(nm,' '));
op = optimset('Display','iter','MaxFunEvals',maxfe,'MaxIter',maxfe);
[xb,Jb] = fminsearch(@obj, x0, op);
pb = setp(xb);
R = struct('pa',pb,'J0',J0,'J',Jb,'x',xb,'names',{nm});
for i=1:6, fprintf('  %-4s %12.5g  ->  %12.5g\n', nm{i}, v0(i), pb.(nm{i})); end
fprintf('\n  maperr %.4f -> %.4f\n', J0, Jb);
% d3 resonance vs BM resonance, the quantity the physics discussion cares about
f1 = sqrt(pb.k1o/pb.m1o)/(2*pi); f5 = sqrt(pb.k5o/pb.m5o)/(2*pi);
fprintf('  BM f1 %.0f Hz | DOF3 f5 %.0f Hz | f5 is %.2f oct below f1 at base\n', ...
        f1, f5, log2(f1/f5));
fprintf('  slope match: k1e-m1e = %.3f, k5e-m5e = %.3f (equal = fixed interval)\n', ...
        pb.k1e-pb.m1e, pb.k5e-pb.m5e);
save('fit_dof3.mat','R');
end
