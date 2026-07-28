C=load('parfit26_c4map_m2b.mat.ckpt.mat');
f=fieldnames(C); if(numel(f)==1 && isstruct(C.(f{1}))), C=C.(f{1}); end
fn=fieldnames(C);
fprintf('checkpoint fields: %s\n', strjoin(fn',', '));
for k={'iter','fev','J','when'}
    if(isfield(C,k{1}))
        v=C.(k{1});
        if(ischar(v)), fprintf('  %-5s = %s\n', k{1}, v);
        else, fprintf('  %-5s = %g\n', k{1}, v); end
    end
end
disp('CKPT_READ_DONE');
