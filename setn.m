function pa = setn(pa, n)
%SETN  Change the spatial grid size, resampling every n-length field.
%
%   pa = setn(pa, n)
%
%   Changing pa.n alone is NOT enough: the parameter struct carries PROFILE
%   ARRAYS defined on the old grid, and leaving them at the old length throws
%   "Arrays have incompatible sizes for this operation" the moment tdm26 builds
%   the cochlea. As of 2026-07-25 those fields are:
%       pa.gampro   cochlear-amplifier gain profile   (tdm26.m ~1808)
%       pa.synpro   synaptic/neural drive profile
%   Both are present in a FRESH modpar26(k) as well as in saved fits, so this
%   affects every configuration.
%
%   Profiles are resampled by linear interpolation over NORMALIZED position, so
%   a non-uniform profile keeps its shape. pa.isv (save locations) is likewise
%   rescaled by fraction, since its entries are indices into 1..n.
%
%   The field list is discovered by LENGTH rather than hardcoded, so a profile
%   added later is handled automatically. That matters because a silently
%   unresampled profile does not error in every code path; it can also just
%   apply the wrong values.

nold = pa.n;
if (n == nold), return; end
xo = linspace(0, 1, nold);
xn = linspace(0, 1, n);

fn = fieldnames(pa);
for i = 1:numel(fn)
    v = pa.(fn{i});
    if (isnumeric(v) && numel(v) == nold && ~strcmp(fn{i},'isv'))
        col = size(v,1) == nold;                 % preserve orientation
        vi = interp1(xo, v(:), xn, 'linear', 'extrap');
        if (col), pa.(fn{i}) = vi(:); else, pa.(fn{i}) = vi(:).'; end
    end
end

if (isfield(pa,'isv') && ~isempty(pa.isv))       % indices -> same fractions
    pa.isv = unique(max(1, min(n, round((pa.isv/nold)*n))), 'stable');
end
pa.n = n;
end
