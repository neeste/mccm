% Targeted test script for specific ABR latency discrepancies

function test_tweaks()
    try
        % Get baseline parameters (1-chamber)
        sav = tdm26(0, 1);
        pa_base = sav.pa;
        close all; % close any figures opened
        
        % Test 1: Baseline
        fprintf('\n--- Baseline ---\n');
        run_tbabr(pa_base);
        
        % Test 2: Decrease ihcex to 4 (extreme)
        pa_test = pa_base;
        pa_test.ihcex = 4;
        fprintf('\n--- Test: ihcex = 4 ---\n');
        run_tbabr(pa_test);

        % Test 3: Increase ihcdr by 10x
        pa_test = pa_base;
        pa_test.ihcdr = pa_base.ihcdr * 10;
        fprintf('\n--- Test: ihcdr = 10x ---\n');
        run_tbabr(pa_test);
        
        % Test 4: ihctc (time constant) = 0.1 ms (half)
        pa_test = pa_base;
        pa_test.ihctc = 0.1e-3;
        fprintf('\n--- Test: ihctc = 0.1ms ---\n');
        run_tbabr(pa_test);
    catch e
        fprintf('FATAL ERROR: %s\n', e.message);
        for k=1:length(e.stack)
            fprintf('  in %s (line %d)\n', e.stack(k).name, e.stack(k).line);
        end
    end
end

function run_tbabr(pa)
    diary('tweak_log.txt');
    try
        tdm26('tbabr', pa);
    catch e
        fprintf('ERROR running tbabr: %s\n', e.message);
    end
    diary off;
    
    fid = fopen('tweak_log.txt', 'r');
    if fid == -1
        fprintf('Could not open tweak_log.txt\n');
        return;
    end
    while ~feof(fid)
        line = fgetl(fid);
        if startsWith(line, ' 0.50  80')
            fprintf('Result -> 0.5 kHz @ 80 dB latency: %s\n', strtrim(line));
        end
        if startsWith(line, ' 4.00  20')
            fprintf('Result -> 4.0 kHz @ 20 dB latency: %s\n', strtrim(line));
        end
    end
    fclose(fid);
    delete('tweak_log.txt');
end
