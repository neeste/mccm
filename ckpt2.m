C=load('parfit26_c3amp.mat.ckpt.mat');
f=fieldnames(C); if(numel(f)==1 && isstruct(C.(f{1}))), C=C.(f{1}); end
fprintf('iter=%g  fevals=%g  J=%.4f  when=%s\n', C.iter, C.fevals, C.J, C.when);
if (C.J >= 1e5)
    fprintf('*** J IS THE 1e6 SENTINEL -- the guard is rejecting everything. KILL IT. ***\n');
else
    fprintf('J is finite and below sentinel -- the objective is live.\n');
end
disp('CKPT2_DONE');
