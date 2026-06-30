clear; clc; close all

regime    = 'b';
pert_type = 'local';
k_mode    = 2;
pert_mag  = 0.05;
t_end     = 0.02;

IC_list = [
    0.2, 0.2;
    0.2, 0.2;
];


p = getParameters(regime);

C1 = [0.0000 0.0000 1.0000];
C2 = [1.0000 0.0000 0.0000];

regime_names = struct('b', 'Bistable', 'o', 'Oscillatory', 'c', 'Coexistence');
regime_label = regime_names.(regime);

folder_name = ['figure_' regime '_' pert_type];
if strcmp(pert_type, 'cosine')
    folder_name = [folder_name '_k' num2str(k_mode)];
end
if ~exist(folder_name, 'dir')
    mkdir(folder_name);
end

initial_guess = [0.5; 0.5];
options_fs = optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-6, 'TolX', 1e-6);
[ss_vals, ~, exitflag] = fsolve(@(y) Model3_ODE_internal(0, y, p), initial_guess, options_fs);

if exitflag <= 0
    [~, y_int] = ode45(@(t,y) Model3_ODE_internal(t, y, p), [0 1000], initial_guess);
    ss_vals = y_int(end, :)';
end

fprintf('\n=== Steady State Found ===\n');
fprintf('  R = %.4f,  rho = %.4f\n\n', ss_vals(1), ss_vals(2));

set(groot, ...
    'defaultFigureColor',          'w', ...
    'defaultAxesColor',            'w', ...
    'defaultAxesFontSize',         12,  ...
    'defaultAxesFontName',         'Times New Roman', ...
    'defaultAxesLineWidth',        0.8, ...
    'defaultAxesBox',              'off', ...
    'defaultAxesTickDir',          'out', ...
    'defaultLineLineWidth',        1.5, ...
    'defaultTextFontName',         'Times New Roman', ...
    'defaultTextFontSize',         12,  ...
    'defaultLegendBox',            'off', ...
    'defaultLegendFontSize',       11);

x  = linspace(0, 1, 500);
t  = linspace(0, t_end, 1000);
L  = x(end) - x(1);
Nx = length(x);
dt = diff(t);

nIC = size(IC_list, 1);

all_sol = cell(nIC, 1);
all_L2  = cell(nIC, 1);

a5_str   = strrep(num2str(p.a5,  '%.4g'), '.', 'p');
tend_str = strrep(num2str(t_end, '%.4g'), '.', 'p');

for icIdx = 1:nIC

    dR_manual   = IC_list(icIdx, 1);
    drho_manual = IC_list(icIdx, 2);

    ic_opts.pert_type   = pert_type;
    ic_opts.pert_mag    = pert_mag;
    ic_opts.k_mode      = k_mode;
    ic_opts.dR_manual   = dR_manual;
    ic_opts.drho_manual = drho_manual;

    fprintf('--- Running PDE for IC %d:  dR = %.4f,  drho = %.4f ---\n', ...
            icIdx, dR_manual, drho_manual);

    sol = pdepe(0, ...
        @(x,t,u,dudx) Model3_PDE(x,t,u,dudx,p), ...
        @(x) Model3_IC(x, ss_vals, ic_opts, Nx), ...
        @Model3_BC, x, t);

    all_sol{icIdx} = sol;

    % L2 norm of du/dt
    L2 = zeros(length(t)-1, 4);
    for i = 1:length(t)-1
        for j = 1:4
            u_diff = (sol(i+1,:,j) - sol(i,:,j)) / dt(i);
            L2(i,j) = sqrt(trapz(x, u_diff.^2));
        end
    end
    L2_range = max(L2) - min(L2);
    L2_range(L2_range == 0) = 1;
    all_L2{icIdx} = (L2 - min(L2)) ./ L2_range;

    file_prefix = sprintf('%s_a5_%s_t_%s_%s_IC%d', ...
                          regime_label, a5_str, tend_str, pert_type, icIdx);

    % Spatial profiles
    fig1 = figure('Name', sprintf('Spatial Profiles IC%d', icIdx), 'Color', 'w');
    subplot(1,2,1);
    plot(x, sol(end,:,1), 'Color', C1, 'LineWidth', 1.5);
    title('Active Rac (R)'); xlabel('x'); ylabel('Concentration');
    applyThesisAxes(gca);
    subplot(1,2,2);
    plot(x, sol(end,:,3), 'Color', C2, 'LineWidth', 1.5);
    title('Active Rho (\rho)'); xlabel('x'); ylabel('Concentration');
    applyThesisAxes(gca);
    sgtitle(sprintf('IC %d  |  (\\DeltaR=%.2g, \\Delta\\rho=%.2g)  |  %s', ...
            icIdx, dR_manual, drho_manual, pert_type), ...
            'FontName','Times New Roman','FontSize',13,'FontWeight','bold');
    exportgraphics(fig1, fullfile(folder_name, [file_prefix '_spatial.png']), ...
                   'Resolution', 300, 'BackgroundColor', 'white');
    close(fig1);

    % L2 convergence
    fig2 = figure('Name', sprintf('L2 IC%d', icIdx), 'Color', 'w'); hold on;
    plot(t(1:end-1), all_L2{icIdx}(:,1), 'Color', C1, 'LineWidth', 1.5, 'DisplayName', 'R');
    plot(t(1:end-1), all_L2{icIdx}(:,3), 'Color', C2, 'LineWidth', 1.5, 'DisplayName', '\rho');
    title('Normalised L^2 Norm of \partialu/\partialt');
    xlabel('Time'); ylabel('$\|.|_{L^2}$', 'Interpreter', 'latex');
    legend('Location', 'northeast'); xlim([0, t_end]);
    applyThesisAxes(gca);
    exportgraphics(fig2, fullfile(folder_name, [file_prefix '_L2.png']), ...
                   'Resolution', 300, 'BackgroundColor', 'white');
    close(fig2);

    % Surfaces
    fig4 = figure('Name', sprintf('Surfaces IC%d', icIdx), 'Color', 'w', ...
                  'Position', [100 100 1400 500]);
    subplot(1,2,1);
    pcolor(x, t, sol(:,:,1)); shading interp; colormap(gca, 'parula'); axis tight;
    cb1 = colorbar; cb1.Label.String = 'Concentration';
    cb1.Label.FontName = 'Times New Roman'; cb1.Label.FontSize = 11;
    title('Active Rac (R)'); xlabel('x'); ylabel('t');
    applyThesisAxes(gca);
    subplot(1,2,2);
    pcolor(x, t, sol(:,:,3)); shading interp; colormap(gca, 'parula'); axis tight;
    cb2 = colorbar; cb2.Label.String = 'Concentration';
    cb2.Label.FontName = 'Times New Roman'; cb2.Label.FontSize = 11;
    title('Active Rho (\rho)'); xlabel('x'); ylabel('t');
    applyThesisAxes(gca);
    sgtitle(sprintf('IC %d  |  (\\DeltaR=%.2g, \\Delta\\rho=%.2g)  |  %s', ...
            icIdx, dR_manual, drho_manual, pert_type), ...
            'FontName','Times New Roman','FontSize',13,'FontWeight','bold');
    exportgraphics(fig4, fullfile(folder_name, [file_prefix '_surface.png']), ...
                   'Resolution', 300, 'BackgroundColor', 'white');
    close(fig4);

    % Mass conservation
    init_Rac  = (1/L) * trapz(x, sol(1,:,1)  + sol(1,:,2));
    init_Rho  = (1/L) * trapz(x, sol(1,:,3)  + sol(1,:,4));
    total_Rac = (1/L) * trapz(x, sol(:,:,1)  + sol(:,:,2), 2);
    total_Rho = (1/L) * trapz(x, sol(:,:,3)  + sol(:,:,4), 2);
    fig3 = figure('Name', sprintf('Mass IC%d', icIdx), 'Color', 'w');
    subplot(1,2,1);
    plot(t, total_Rac, 'Color', C1, 'LineWidth', 1.5); hold on;
    yline(init_Rac, 'k--', 'LineWidth', 1, 'Label', 'Initial');
    title('Total Rac'); xlabel('Time'); ylabel('Mean concentration');
    ylim([init_Rac-0.5, init_Rac+0.5]); applyThesisAxes(gca);
    subplot(1,2,2);
    plot(t, total_Rho, 'Color', C2, 'LineWidth', 1.5); hold on;
    yline(init_Rho, 'k--', 'LineWidth', 1, 'Label', 'Initial');
    title('Total Rho'); xlabel('Time'); ylabel('Mean concentration');
    ylim([init_Rho-0.5, init_Rho+0.5]); applyThesisAxes(gca);
    sgtitle('Mass Conservation', 'FontName','Times New Roman','FontSize',13,'FontWeight','bold');
    exportgraphics(fig3, fullfile(folder_name, [file_prefix '_mass.png']), ...
                   'Resolution', 300, 'BackgroundColor', 'white');
    close(fig3);

    % IC profile
    IC_all = zeros(Nx, 4);
    for ix = 1:Nx
        IC_all(ix,:) = Model3_IC(x(ix), ss_vals, ic_opts, Nx)';
    end
    fig5 = figure('Name', sprintf('IC Profile %d', icIdx), 'Color', 'w');
    plot(x, IC_all(:,1), 'Color', C1, 'LineWidth', 1.5, 'DisplayName', 'R'); hold on;
    plot(x, IC_all(:,3), 'Color', C2, 'LineWidth', 1.5, 'DisplayName', '\rho');
    legend('Location', 'best'); xlabel('x'); ylabel('Concentration');
    title(sprintf('IC %d — %s  (\\DeltaR=%.2g, \\Delta\\rho=%.2g)', ...
          icIdx, pert_type, dR_manual, drho_manual));
    applyThesisAxes(gca);
    exportgraphics(fig5, fullfile(folder_name, [file_prefix '_IC.png']), ...
                   'Resolution', 300, 'BackgroundColor', 'white');
    close(fig5);
end

nCols      = 3;
nRowsEach  = 2;          % ICs per image
panelW     = 420;
panelH     = 280;
figW       = nCols * panelW;
figH       = nRowsEach * panelH;

panel_labels = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
nImages      = ceil(nIC / nRowsEach);   % = 2 for 4 ICs

for imgIdx = 1:nImages

    ic_start = (imgIdx - 1) * nRowsEach + 1;
    ic_end   = min(imgIdx * nRowsEach, nIC);
    nRowsThis = ic_end - ic_start + 1;

    figComb = figure('Name', sprintf('Combined Panel %d', imgIdx), 'Color', 'w', ...
                     'Position', [50 50 figW nRowsThis*panelH]);

    labelIdx = (imgIdx - 1) * nRowsEach * nCols + 1;   % A=1, G=7, etc.
    localRow = 0;

    for icIdx = ic_start:ic_end
        localRow = localRow + 1;
        sol_i = all_sol{icIdx};
        L2_i  = all_L2{icIdx};

        % --- Column 1: Active Rac surface ---
        ax1 = subplot(nRowsThis, nCols, (localRow-1)*nCols + 1);
        pcolor(x, t, sol_i(:,:,1)); shading interp; colormap(ax1, 'parula'); axis tight;
        cb = colorbar;
        cb.Label.String   = 'Conc.';
        cb.Label.FontName = 'Times New Roman';
        cb.Label.FontSize = 10;
        xlabel('x'); ylabel('t');
        title(ax1, ['\bf' panel_labels(labelIdx)], ...
              'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'bold');
        applyThesisAxes(ax1);
        labelIdx = labelIdx + 1;

        % --- Column 2: Active Rho surface ---
        ax2 = subplot(nRowsThis, nCols, (localRow-1)*nCols + 2);
        pcolor(x, t, sol_i(:,:,3)); shading interp; colormap(ax2, 'parula'); axis tight;
        cb = colorbar;
        cb.Label.String   = 'Conc.';
        cb.Label.FontName = 'Times New Roman';
        cb.Label.FontSize = 10;
        xlabel('x'); ylabel('t');
        title(ax2, ['\bf' panel_labels(labelIdx)], ...
              'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'bold');
        applyThesisAxes(ax2);
        labelIdx = labelIdx + 1;

        % --- Column 3: L2 norm ---
        ax3 = subplot(nRowsThis, nCols, (localRow-1)*nCols + 3);
        hold(ax3, 'on');
        plot(ax3, t(1:end-1), L2_i(:,1), 'Color', C1, 'LineWidth', 1.5, 'DisplayName', 'R');
        plot(ax3, t(1:end-1), L2_i(:,3), 'Color', C2, 'LineWidth', 1.5, 'DisplayName', '\rho');
        xlabel(ax3, 'Time');
        ylabel(ax3, '$\|\partial u / \partial t\|_{L^2}$', 'Interpreter', 'latex');
        title(ax3, ['\bf' panel_labels(labelIdx)], ...
              'FontName', 'Times New Roman', 'FontSize', 14, 'FontWeight', 'bold');
        legend(ax3, 'Location', 'northeast');
        xlim(ax3, [0, t_end]);
        applyThesisAxes(ax3);
        labelIdx = labelIdx + 1;
    end

    combined_file = fullfile(folder_name, ...
        sprintf('%s_a5_%s_t_%s_%s_combined_part%d.png', ...
                regime_label, a5_str, tend_str, pert_type, imgIdx));
    exportgraphics(figComb, combined_file, 'Resolution', 300, 'BackgroundColor', 'white');
    close(figComb);
    fprintf('Combined panel %d saved to: %s\n', imgIdx, combined_file);
end

% Reset groot defaults
set(groot, ...
    'defaultFigureColor',      'default', ...
    'defaultAxesColor',        'default', ...
    'defaultAxesFontSize',     'default', ...
    'defaultAxesFontName',     'default', ...
    'defaultAxesLineWidth',    'default', ...
    'defaultAxesBox',          'default', ...
    'defaultAxesTickDir',      'default', ...
    'defaultLineLineWidth',    'default', ...
    'defaultTextFontName',     'default', ...
    'defaultTextFontSize',     'default', ...
    'defaultLegendBox',        'default', ...
    'defaultLegendFontSize',   'default');

% =========================================================
%  LOCAL FUNCTIONS
% =========================================================

function applyThesisAxes(ax)
    set(ax, ...
        'Color',     'w', ...
        'Box',       'off', ...
        'TickDir',   'out', ...
        'FontSize',  11, ...
        'FontName',  'Times New Roman', ...
        'LineWidth', 0.8, ...
        'GridColor', [0.85 0.85 0.85], ...
        'GridAlpha', 1.0);
    grid(ax, 'on');
end

function p = getParameters(regime)
    p.n     = 1;
    p.gamma = 2000;
    p.dv    = 100;
    switch lower(regime)
        case 'b'
            p.a1 = 0.45; p.k1 = 1.0;
            p.a2 = 6.0;  p.k2 = 0.02;
            p.a3 = 5.5;  p.k3 = 0.30;
            p.a4 = 1.6;  p.k4 = 1.5;
            p.k5 = 0.75; p.a5 = 1; p.a6 = 1.0;
        case 'o'
            p.a1 = 0.25; p.k1 = 2.5;
            p.a2 = 5.5;  p.k2 = 0.075;
            p.a3 = 6.0;  p.k3 = 0.30;
            p.a4 = 0.1;  p.k4 = 1.5;
            p.k5 = 0.55; p.a5 = 0.05; p.a6 = 1.0;
        case 'c'
            p.a1 = 0.50; p.k1 = 3.0;
            p.a2 = 7.0;  p.k2 = 0.075;
            p.a3 = 4.7;  p.k3 = 0.25;
            p.a4 = 0.2;  p.k4 = 2.0;
            p.k5 = 0.50; p.a5 = 12;   p.a6 = 0.6;
        otherwise
            error('Use regime = b, o, or c');
    end
end

function dy = Model3_ODE_internal(~, y, p)
    p.gamma = 1;
    R = y(1); rho = y(2);
    dy = [
        (p.a1*(1-R)^p.n)  /(p.k1^p.n+(1-R)^p.n)   + ...
        (p.a2*(1-R)^p.n*R)/(p.k2^p.n+(1-R)^p.n)   - ...
        (p.a3*rho*R^p.n)  /(p.k3^p.n+R^p.n)        - R;

        (p.a4*(1-rho)^p.n)  /(p.k4^p.n+(1-rho)^p.n) + ...
        (p.a5*R*(1-rho)^p.n)/(p.k5^p.n+(1-rho)^p.n) - p.a6*rho
    ];
end

function [c,f,s] = Model3_PDE(~,~,u,dudx,p)
    c = [1;1;1;1];
    f = [1; p.dv; 1; p.dv] .* dudx;
    R   = u(1); Ri   = u(2);
    rho = u(3); rhoi = u(4);
    TermR = (p.a1*Ri^p.n)  /(p.k1^p.n+Ri^p.n)   + ...
            (p.a2*Ri^p.n*R)/(p.k2^p.n+Ri^p.n)   - ...
            (p.a3*rho*R^p.n)/(p.k3^p.n+R^p.n)   - R;
    TermRho = (p.a4*rhoi^p.n)  /(p.k4^p.n+rhoi^p.n) + ...
              (p.a5*R*rhoi^p.n)/(p.k5^p.n+rhoi^p.n) - p.a6*rho;
    s = [p.gamma*TermR; -p.gamma*TermR; p.gamma*TermRho; -p.gamma*TermRho];
end

function u0 = Model3_IC(x, ss_vals, ic_opts, ~)
    R_ss      = ss_vals(1);
    rho_ss    = ss_vals(2);
    dR_base   = ic_opts.dR_manual;
    drho_base = ic_opts.drho_manual;
    switch lower(ic_opts.pert_type)
        case 'manual'
            R   = R_ss   + dR_base;
            rho = rho_ss + drho_base;
        case 'random'
            dR   = ic_opts.pert_mag * (rand() - 0.5);
            drho = ic_opts.pert_mag * (rand() - 0.5);
            R   = R_ss   + dR_base + dR;
            rho = rho_ss + drho_base + drho;
        case 'cosine'
            R   = R_ss   + dR_base   + ic_opts.pert_mag * cos(ic_opts.k_mode * pi * x);
            rho = rho_ss + drho_base - ic_opts.pert_mag * cos(ic_opts.k_mode * pi * x);
        case 'local'
            x_c   = 0.5;
            sigma = 0.05;
            bump  = exp(-((x - x_c).^2) / (2*sigma^2));
            R   = R_ss   + dR_base   + ic_opts.pert_mag * bump;
            rho = rho_ss + drho_base - ic_opts.pert_mag * bump;
        otherwise
            error('pert_type must be ''manual'', ''random'', ''cosine'', or ''local''.');
    end
    R   = max(0.001, min(0.999, R));
    rho = max(0.001, min(0.999, rho));
    u0  = [R; 1-R; rho; 1-rho];
end

function [pl,ql,pr,qr] = Model3_BC(~,~,~,~,~)
    pl = zeros(4,1);
    ql = ones(4,1);
    pr = zeros(4,1);
    qr = ones(4,1);
end