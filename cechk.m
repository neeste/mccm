% CECHK -- acceptance tests for coupeig on [d;v;vohc].
% 1. OFF UNCHANGED: nv=0 must reproduce the committed spectrum exactly.
% 2. THE POLE NOW REGISTERS: maxRe/maxRe_osc must MOVE with fc. Before this they
%    were identical at fc = Inf/8000/4000/1000, which is what "blind" looked like.
% 3. THE RC EIGENVALUES ARE THERE AND ARE RIGHT: an isolated one-pole block
%    contributes eigenvalues near -1/tau. Their presence at the predicted place
%    is the direct evidence the new rows are wired correctly, not just that
%    something changed.
% 4. OPERATOR SIZE is 2*nsv + n when on, 2*nsv when off.
PROJ = '/Users/neely/Library/CloudStorage/OneDrive-FatherFlanagan''sBoysHome/STN/mod/mccm';
addpath(PROJ); cd(PROJ);
L=load('fit_nch3_surface.mat'); pa=L.R.pa;
if (isfield(pa,'hbmode')), pa=rmfield(pa,'hbmode'); end
fprintf('\n=== 1+2+3: spectrum vs corner frequency ===\n');
fprintf('  %9s %8s %11s %12s %14s %10s\n','fc (Hz)','N','maxRe','maxRe_osc','-1/tau','nearest lam');
for fc=[Inf 8000 4000 1000 500]
    p=pa; if (isfinite(fc)), p.ohctau=1/(2*pi*fc); end
    try
        evalc('C=tdm26(''coupeig'',struct(''pa'',p));');
        N=numel(C.lam);
        if (isfinite(fc))
            invtau=-1/p.ohctau;
            [~,j]=min(abs(C.lam-invtau));
            fprintf('  %9.4g %8d %+11.3f %+12.3f %14.1f %10.1f\n', ...
                    fc, N, C.maxRe, C.maxRe_osc, invtau, real(C.lam(j)));
        else
            fprintf('  %9.4g %8d %+11.3f %+12.3f %14s %10s\n', fc, N, C.maxRe, C.maxRe_osc, '-', '-');
        end
    catch e
        fprintf('  %9.4g FAILED: %s\n', fc, e.message(1:min(60,end)));
    end
end
fprintf('\n  OFF baseline expected: maxRe +22.07, maxRe_osc -180.23\n');
fprintf('  If maxRe/maxRe_osc are now IDENTICAL across fc, the extension did not take.\n');
fprintf('\n=== 4: operator size ===\n');
p=pa; p.ohctau=1/(2*pi*1000);
d=resolve_dof_probe(p); nsv=p.n*d+max(p.nmev,1);
fprintf('  n %d dof %d nsv %d -> expect N off %d, N on %d\n', p.n, d, nsv, 2*nsv, 2*nsv+p.n);
disp('CECHK_DONE');
function d=resolve_dof_probe(pa)
d=2; if (pa.m>=4 || (isfield(pa,'d3int')&&pa.d3int)), d=3; end
end
