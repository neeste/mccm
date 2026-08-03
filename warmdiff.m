% WARMDIFF -- which fields of a warm start do NOT survive into jointobj?
%
% parfit26 REPORTS its start from pa0 (the warm start) but EVALUATES
% setpar_l(base,pv) with base=modpar26(nch). Any tuned field outside the 30
% parnames entries + chsz (+ gpgrade/hbsc) is therefore replaced by a modpar26
% default on every evaluation, while the printed start line still shows the good
% model. Measured consequence on the running nch=3 fit: the start line says
% maperr 147.1 and the very first jointobj evaluation scores maperr ~934.
%
% This is invisible whenever the warm start IS a modpar26 baseline -- which is
% why the nch=1 runs were unaffected: sweep1b fell back to modpar26(1), so
% base+pv reproduced pa0 exactly.
H = parfit26('handles');
for f = {'parfit26_nch3.mat'}
    L = load(f{1}); pa0 = L.R.pa;
    base = modpar26(3);
    pv = H.getpar(pa0);
    par = H.setpar(base, pv, []);      % exactly what jointobj evaluates
    fprintf('\n== %s: fields lost between pa0 and setpar_l(modpar26(3),pv) ==\n', f{1});
    nm = union(fieldnames(pa0), fieldnames(base));
    nlost = 0;
    for i = 1:numel(nm)
        k = nm{i};
        a = []; b = [];
        if (isfield(pa0,k)), a = pa0.(k); end
        if (isfield(par,k)), b = par.(k); end
        if (~isnumeric(a) || ~isnumeric(b)), continue; end
        if (isempty(a) && isempty(b)), continue; end
        same = isequal(size(a),size(b)) && (isempty(a) || max(abs(double(a(:))-double(b(:)))) < 1e-12);
        if (~same)
            nlost = nlost + 1;
            if (isscalar(a) && isscalar(b))
                fprintf('   %-12s  pa0 %-14.6g  jointobj %-14.6g\n', k, a, b);
            else
                fprintf('   %-12s  pa0 [%s] range [%.4g %.4g]   jointobj [%s] range [%.4g %.4g]\n', ...
                    k, mat2str(size(a)), min(a(:)), max(a(:)), mat2str(size(b)), min(b(:)), max(b(:)));
            end
        end
    end
    if (nlost==0), fprintf('   (none -- warm start survives intact)\n'); end
    fprintf('   %d field(s) differ\n', nlost);
end
disp('WARMDIFF_DONE');
