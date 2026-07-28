% Which pa fields are n-length arrays? Those must be regenerated when n changes.
L=load('refit_c1broad.mat'); pa=L.R.pa;
fn=fieldnames(pa); n=pa.n;
fprintf('pa.n = %d\n\nn-length or other non-scalar fields:\n', n);
for i=1:numel(fn)
    v=pa.(fn{i});
    if (isnumeric(v) && numel(v)>1)
        tag=''; if (numel(v)==n), tag='  <== n-LENGTH'; end
        fprintf('  %-12s numel=%-6d %s\n', fn{i}, numel(v), tag);
    end
end
% same for a fresh modpar26(1)
p2=modpar26(1); fn2=fieldnames(p2);
fprintf('\nfresh modpar26(1) non-scalar fields:\n');
for i=1:numel(fn2)
    v=p2.(fn2{i});
    if (isnumeric(v) && numel(v)>1)
        tag=''; if (numel(v)==p2.n), tag='  <== n-LENGTH'; end
        fprintf('  %-12s numel=%-6d %s\n', fn2{i}, numel(v), tag);
    end
end
disp('FIND_NARR_DONE');
