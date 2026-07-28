% Generate baseline presentation figures without the slow fwdmsk
try
    fprintf('Running m=1 tbabr...\n');
    tdm26('tbabr', 1, 4, 2);   % args 3&4 nonzero -> generate the latency figure
    saveas(figure(1), 'm1_tbabr.png');

    fprintf('Running m=3 tbabr...\n');
    tdm26('tbabr', 3, 4, 2);   % args 3&4 nonzero -> generate the latency figure
    saveas(figure(1), 'm3_tbabr.png');
    
    fprintf('Running 3-chamber features...\n');
    tdm26(0, 3);
    figs = findobj('Type','figure');
    for i=1:length(figs)
        saveas(figs(i), sprintf('m3_tdm_feat%d.png', figs(i).Number));
    end
    
    fdm26(3);
    figs = findobj('Type','figure');
    for i=1:length(figs)
        saveas(figs(i), sprintf('m3_fdm_feat%d.png', figs(i).Number));
    end
catch e
    fprintf('Error: %s\n', e.message);
end
exit;
