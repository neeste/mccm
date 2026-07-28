function plot_figs_2008(nch)
    if nargin < 1, nch = 3; end

    % Initialize parameters
    pa = tdm26(0, nch); 
    pa = pa.pa;
    pa.gam = 1; % Active cochlea
    
    dx = pa.xl / (pa.n - 1);
    x = transpose(linspace(0, pa.xl, pa.n));
    x_tap = x .* (1 + (pa.xtap * x) .^ pa.xtex);
    q = x_tap .^ 2;
    bw = pa.bwo * exp(pa.bwe * x_tap + pa.bwq * q);
    ac = pa.aco * exp(pa.ace * x_tap + pa.acq * q);
    
    chsz = pa.chsz;
    chsz = chsz * (2 / sum(chsz));
    A1 = ac * chsz(1);
    A3 = ac * chsz(end); % End is the tympanic scala
    
    %% FIGURE 2 (Single frequency wave decomposition)
    f_eval = 1000; % 1 kHz
    s = 2i * pi * f_eval;
    
    Zs = s * pa.rho * (1 ./ A1 + 1 ./ A3); % Effective Series Impedance
    
    % Smooth Model
    if isfield(pa, 'rough_amp'), pa = rmfield(pa, 'rough_amp'); end
    [~, Yd_smooth, Pd_smooth, ~, ~, ~, ~] = fdmod23(pa, f_eval);
    Yb_smooth = Yd_smooth .* bw;
    
    % Rough Model
    pa.rough_amp = 1e-5;
    [~, Yd_rough, Pd_rough, ~, ~, ~, ~] = fdmod23(pa, f_eval);
    Yb_rough = Yd_rough .* bw;
    
    % Wave decomposition
    [Pp_smooth, Pm_smooth] = extract_waves(Pd_smooth, Yb_smooth, Zs, dx);
    [Pp_rough, Pm_rough]   = extract_waves(Pd_rough, Yb_rough, Zs, dx);
    
    % Plot Figure 2
    figure(201); clf;
    
    % Smooth Magnitude
    subplot(2, 2, 1);
    plot(x, 20*log10(abs(Pp_smooth)), 'k', 'LineWidth', 1.5); hold on;
    plot(x, 20*log10(abs(Pm_smooth)), 'k--', 'LineWidth', 1.5);
    title('Smooth BM'); ylabel('Level (dB)'); xlim([0 35]); ylim([-80 60]); grid on;
    
    % Rough Magnitude
    subplot(2, 2, 2);
    plot(x, 20*log10(abs(Pp_rough)), 'k', 'LineWidth', 1.5); hold on;
    plot(x, 20*log10(abs(Pm_rough)), 'k--', 'LineWidth', 1.5);
    title('Rough BM'); xlim([0 35]); ylim([-80 60]); grid on;
    
    % Smooth Phase
    subplot(2, 2, 3);
    plot(x, unwrap(angle(Pp_smooth))/(2*pi), 'k', 'LineWidth', 1.5); hold on;
    plot(x, unwrap(angle(Pm_smooth))/(2*pi), 'k--', 'LineWidth', 1.5);
    xlabel('Distance from stapes (mm)'); ylabel('Phase (cyc)'); xlim([0 35]); grid on;
    
    % Rough Phase
    subplot(2, 2, 4);
    plot(x, unwrap(angle(Pp_rough))/(2*pi), 'k', 'LineWidth', 1.5); hold on;
    plot(x, unwrap(angle(Pm_rough))/(2*pi), 'k--', 'LineWidth', 1.5);
    xlabel('Distance from stapes (mm)'); xlim([0 35]); grid on;
    
    saveas(gcf, sprintf('/Users/neely/.gemini/antigravity/brain/7ec06545-a0d3-4e3d-8741-3d8d1db66f2f/fig2_m%d.png', nch));
    
    %% FIGURE 3 (Peak Pressure vs Frequency, Exact and Approximate)
    flst = 500 * 2.^linspace(-2, 5.3, 100); % Log spacing from ~125 Hz to 20 kHz
    
    pk_ex_sm = zeros(size(flst));
    pk_ex_ro = zeros(size(flst));
    pk_ap_sm = zeros(size(flst));
    pk_ap_ro = zeros(size(flst));
    
    fprintf('Running frequency sweep for Figure 3...\n');
    for k = 1:length(flst)
        f = flst(k);
        s = 2i * pi * f;
        Zs = s * pa.rho * (1 ./ A1 + 1 ./ A3);
        
        % Smooth
        if isfield(pa, 'rough_amp'), pa = rmfield(pa, 'rough_amp'); end
        [~, Yd, Pd, ~, ~, ~, ~] = fdmod23(pa, f);
        Yb = Yd .* bw;
        pk_ex_sm(k) = 20*log10(max(abs(Pd)));
        Pd_approx = wkb_approx(Pd(1), Yb, Zs, dx);
        pk_ap_sm(k) = 20*log10(max(abs(Pd_approx)));
        
        % Rough
        pa.rough_amp = 1e-5;
        [~, Yd, Pd, ~, ~, ~, ~] = fdmod23(pa, f);
        Yb = Yd .* bw;
        pk_ex_ro(k) = 20*log10(max(abs(Pd)));
        Pd_approx = wkb_approx(Pd(1), Yb, Zs, dx);
        pk_ap_ro(k) = 20*log10(max(abs(Pd_approx)));
    end
    
    figure(202); clf;
    semilogx(flst/1000, pk_ex_sm, 'k--', 'LineWidth', 2); hold on;
    semilogx(flst/1000, pk_ex_ro, 'k', 'LineWidth', 2);
    semilogx(flst/1000, pk_ap_sm, 'r--', 'LineWidth', 1);
    semilogx(flst/1000, pk_ap_ro, 'r', 'LineWidth', 1);
    
    title(sprintf('Figure 3 (m=%d)', nch));
    xlabel('Frequency (kHz)');
    ylabel('Peak pressure (dB re stapes)');
    xlim([0.5 20]);
    xticks([1 2 3 5 10 20]);
    xticklabels({'1','2','3','5','10','20'});
    grid on;
    legend('exact-smooth', 'exact-rough', 'approx-smooth', 'approx-rough', 'Location', 'SouthWest');
    
    saveas(gcf, sprintf('/Users/neely/.gemini/antigravity/brain/7ec06545-a0d3-4e3d-8741-3d8d1db66f2f/fig3_m%d.png', nch));
end

function [Pp, Pm] = extract_waves(Pd, Yb, Zs, dx)
    z0 = sqrt(Zs ./ Yb);
    % Filter out NaNs to calculate gradient
    valid = ~isnan(Pd);
    Pd_clean = Pd; Pd_clean(~valid) = 0;
    
    % Compute volume velocity U = -(1/Zs) * d(Pd)/dx
    dPdx = gradient(Pd_clean, dx);
    U = -(1 ./ Zs) .* dPdx;
    
    Pp = 0.5 * (Pd_clean + z0 .* U);
    Pm = 0.5 * (Pd_clean - z0 .* U);
    
    Pp(~valid) = NaN;
    Pm(~valid) = NaN;
    
    % Normalize to stapes pressure
    stapes_p = Pp(1);
    if stapes_p == 0 || isnan(stapes_p), stapes_p = 1; end
    Pp = Pp / stapes_p;
    Pm = Pm / stapes_p;
end

function Pd_approx = wkb_approx(P0, Yb, Zs, dx)
    kappa = sqrt(-Zs .* Yb);
    z0 = sqrt(Zs ./ Yb);
    
    eps = 0.5 * gradient(log(z0), dx);
    
    % Integrate Eq 13 for P+
    int_P = cumtrapz(eps - kappa) * dx;
    Pp_approx = P0 * exp(int_P);
    
    % Integrate Eq 16 for P-
    int_Pm = cumtrapz(kappa + eps) * dx;
    
    % The integral from x to L (backward integral)
    integrand = eps .* exp(-2 * cumtrapz(kappa) * dx);
    
    % Calculate backward integral
    n = length(kappa);
    back_int = zeros(n, 1);
    for i = 1:n
        back_int(i) = trapz(integrand(i:end)) * dx;
    end
    
    Pm_approx = P0 * exp(int_Pm) .* back_int;
    
    Pd_approx = Pp_approx + Pm_approx;
    
    % Normalize
    stapes_p = Pp_approx(1);
    if stapes_p == 0 || isnan(stapes_p), stapes_p = 1; end
    Pd_approx = Pd_approx / stapes_p;
end
