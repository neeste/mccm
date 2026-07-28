try
    fprintf('Running fdm26(3)...\n');
    fdm26(3);
    figs = findobj('Type', 'figure');
    for i = 1:length(figs)
        saveas(figs(i), sprintf('fdm26_fig%d.png', figs(i).Number));
    end
    close all;

    fprintf('Running tdm26(0,3)...\n');
    tdm26(0,3);
    figs = findobj('Type', 'figure');
    for i = 1:length(figs)
        saveas(figs(i), sprintf('tdm26_fig%d.png', figs(i).Number));
    end
    close all;

    fprintf('VERIFICATION COMPLETED SUCCESSFULLY\n');
catch e
    fprintf('ERROR OCCURRED: %s\n', e.message);
    for k = 1:length(e.stack)
        fprintf('File: %s, Name: %s, Line: %d\n', e.stack(k).file, e.stack(k).name, e.stack(k).line);
    end
end
exit;
