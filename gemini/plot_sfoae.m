function plot_sfoae(nch)
    if nargin < 1, nch = 1; end
    
    sav = tdm26(0, nch);
    pa = sav.pa;
    pa.gam = 1; % Active cochlea
    
    % High resolution frequency list for SFOAE fine structure
    % Linear frequency spacing is better for visualizing constant-delay ripples
    f1 = 500;
    f2 = 4000;
    flst = linspace(f1, f2, 2000); 
    
    [~,~,~,~,Dh,Zc,Gf] = fdmod23(pa, flst);
    Ze = midear_imp(pa, flst, Zc);
    
    % Convert to dB for Ze
    Ze_dB = 20*log10(abs(Ze));
    
    % Calculate 'reflectance' or pressure ripple proxy
    % A simple detrending can highlight the fine structure
    window_size = 50;
    Ze_smooth = movmean(Ze_dB, window_size);
    Ze_ripple = Ze_dB - Ze_smooth;
    
    figure(1); clf;
    
    subplot(2,1,1);
    plot(flst/1000, Ze_dB, 'b', 'LineWidth', 1);
    hold on;
    plot(flst/1000, Ze_smooth, 'r--', 'LineWidth', 2);
    title(sprintf('Eardrum Impedance Magnitude (m = %d)', nch));
    ylabel('|Ze| (dB)');
    xlim([f1/1000, f2/1000]);
    grid on;
    
    subplot(2,1,2);
    plot(flst/1000, Ze_ripple, 'k', 'LineWidth', 1);
    title('Ze Fine Structure (Detrended)');
    xlabel('Frequency (kHz)');
    ylabel('\Delta |Ze| (dB)');
    xlim([f1/1000, f2/1000]);
    % Set y-limits tight to see if it's perfectly flat or rippled
    ylim([-0.5 0.5]);
    grid on;
    
    fname = sprintf('/Users/neely/.gemini/antigravity/brain/7ec06545-a0d3-4e3d-8741-3d8d1db66f2f/sfoae_m%%d.png', nch);
    saveas(gcf, fname);
    fprintf('Saved %s\n', fname);
end
