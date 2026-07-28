for m=[2 3 4]
  pa=modpar26(m);
  fprintf('m=%d  n=%d  numel(isv)=%d  isv=%s\n  x/L=%s\n', ...
     m, pa.n, numel(pa.isv), mat2str(pa.isv(:)'), num2str(pa.isv(:)'/pa.n,'%6.3f'));
end
