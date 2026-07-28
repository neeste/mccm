try
    fprintf('Starting optimization for m=1...\n');
    res1 = tdm26('optimize', 1);
    save('opt_m1_result.mat', 'res1');
    fprintf('Done m=1. Starting m=3...\n');
    res3 = tdm26('optimize', 3);
    save('opt_m3_result.mat', 'res3');
    fprintf('Optimization complete!\n');
catch e
    fprintf('Error: %s\n', e.message);
    for i=1:length(e.stack)
        fprintf('  In %s (line %d)\n', e.stack(i).name, e.stack(i).line);
    end
end
exit;
