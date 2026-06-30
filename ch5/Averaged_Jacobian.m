clear; clc; close all
format long g

%% Parameters
scenario = 'S3';   % 'S1', 'S2', or 'S3'

switch scenario
    case 'S1'
        p.a1 = 0.45; p.k1 = 1.0;
        p.a2 = 6.0;  p.k2 = 0.02;
        p.a3 = 5.5;  p.k3 = 0.30;
        p.a4 = 1.6;  p.k4 = 1.5;
        p.k5 = 0.75; p.a6 = 1.0;
        p.n  = 1;    gamma = 1;

    case 'S2'
        p.a1 = 0.25; p.k1 = 2.5;
        p.a2 = 5.5;  p.k2 = 0.075;
        p.a3 = 6.0;  p.k3 = 0.30;
        p.a4 = 0.1;  p.k4 = 1.5;
        p.k5 = 0.55; p.a6 = 1.0;
        p.n  = 1;    gamma = 1;

    case 'S3'
        p.a1 = 0.50; p.k1 = 3.0;
        p.a2 = 7.0;  p.k2 = 0.075;
        p.a3 = 4.7;  p.k3 = 0.25;
        p.a4 = 0.2;  p.k4 = 2.0;
        p.k5 = 0.50; p.a6 = 0.6;
        p.n  = 1;    gamma = 1;

    otherwise
        error('Invalid scenario. Choose S1, S2, or S3.')
end

p.a5  = 4.3;
gamma = 140;

fprintf('Scenario: %s\n', scenario)

%% Output folder
save_folder = fullfile(pwd, scenario);
if ~exist(save_folder, 'dir'), mkdir(save_folder); end

param_tag = strrep(sprintf('%s_a5_%g_gamma_%g', scenario, p.a5, gamma), '.', 'p');

dv_range = 0:0.05:25;
n_wave   = 0.01:0.001:4;
d_values_to_plot = [1, 3.4, 5];
tol = 1e-7;

%% Steady state
options = optimoptions('fsolve', 'Display', 'none');
[ss, ~, exitflag] = fsolve(@(x) rho_kinetics_reduced(x, p), [0.9, 0.7], options);
if exitflag <= 0
    warning('fsolve did not converge.')
end

R_ss    = ss(1);
rho_ss  = ss(2);
Ri_ss   = 1 - R_ss;
rhoi_ss = 1 - rho_ss;

fprintf('Steady state: R=%.6f, Ri=%.6f, rho=%.6f, rhoi=%.6f\n', ...
    R_ss, Ri_ss, rho_ss, rhoi_ss)

%% Stability check
Jred       = get_reduced_jacobian_numeric([R_ss; rho_ss], p, gamma);
lambda_red = eig(Jred);
fprintf('Eigenvalues of reaction Jacobian:\n'); disp(lambda_red)

if max(real(lambda_red)) > tol
    fprintf('Homogeneous steady state is UNSTABLE — using averaged Jacobian.\n\n')

    delta = min([1e-1, 0.25*Ri_ss, 0.25*rhoi_ss]);
    if delta <= 0, delta = 1e-3; end
    x0_prep = [R_ss+delta; Ri_ss-delta; rho_ss+delta; rhoi_ss-delta];

    opts = odeset('RelTol',1e-8, 'AbsTol',1e-10, 'NonNegative',1:4);

    [~, x_prep] = ode15s(@(t,x) Nondim_full(t,x,p,gamma), 0:0.1:500, x0_prep, opts);
    [t_full, x_full] = ode15s(@(t,x) Nondim_full(t,x,p,gamma), 0:0.001:250, x_prep(end,:).', opts);

    [x_limit, t_limit, T_final, deviation] = compute_limit_cycle_4var(t_full, x_full);
    if isempty(x_limit), error('Limit cycle not detected.'); end

    fprintf('Limit cycle: T=%.6f, deviation=%.2e\n\n', T_final, deviation)
    Jac_used = average_jacobian_on_cycle(t_limit, x_limit, p, gamma);

else
    fprintf('Homogeneous steady state is STABLE — using steady-state LSA.\n\n')
    Jac_used = get_4x4_jacobian_full(R_ss, Ri_ss, rho_ss, rhoi_ss, p, gamma);
end

%% Eigenvalue scan over diffusion ratio d and wavenumber n
max_eig_vs_d = zeros(length(dv_range), 1);
all_max_eigs = zeros(length(dv_range), length(n_wave));

for i = 1:length(dv_range)
    dv = dv_range(i);
    for j = 1:length(n_wave)
        k2 = (dv > 0) * (n_wave(j)^2) * (pi^2);
        D  = diag([1, dv, 1, dv]);
        M  = Jac_used - k2*D;
        all_max_eigs(i,j) = max(real(eig(M)));
    end
    max_eig_vs_d(i) = max(all_max_eigs(i,:));
end

%% Critical diffusion value
sign_d  = sign(max_eig_vs_d);
cross_d = find(diff(sign_d) > 0, 1);

if ~isempty(cross_d)
    i1 = cross_d; i2 = cross_d + 1;
    d_crit = dv_range(i1) - max_eig_vs_d(i1) * ...
             (dv_range(i2)-dv_range(i1)) / (max_eig_vs_d(i2)-max_eig_vs_d(i1));
    fprintf('Critical diffusion ratio d_c = %.6f\n', d_crit)
else
    fprintf('No critical crossing detected.\n')
end

%% Figure 1 — d-scan
fig1 = figure(1);
set(fig1, 'Color','w', 'Units','centimeters', 'Position',[2 2 16 12])

plot(dv_range, max_eig_vs_d, 'b-', 'LineWidth',2); hold on
yline(0, '--k', 'LineWidth',1.2)
if ~isempty(cross_d)
    xline(d_crit, '--r', 'LineWidth',1.5)
    text(d_crit+0.3, max(max_eig_vs_d)*0.9, ...
        sprintf('$\\mathbf{d_c = %.2f}$', d_crit), ...
        'Interpreter','latex', 'FontSize',14, 'Color',[1 0 0], ...
        'FontName','Times New Roman', 'HorizontalAlignment','left')
end
xlabel('Diffusion ratio, $d$',         'Interpreter','latex', 'FontSize',14)
ylabel('$\max\,\mathrm{Re}(\lambda)$', 'Interpreter','latex', 'FontSize',14)
grid on; box on
set(gca, 'FontSize',14, 'FontName','Times New Roman', ...
    'TickLabelInterpreter','latex', 'LineWidth',1.0, 'GridAlpha',0.3)

exportgraphics(fig1, fullfile(save_folder, ...
    sprintf('%s_d_scan_max_eigenvalue.jpg', param_tag)), 'Resolution',400)

%% Figure 2 — Dispersion relation
fig2 = figure(2);
set(fig2, 'Color','w', 'Units','centimeters', 'Position',[2 2 18 13])
ax = axes('Parent', fig2);
hold(ax, 'on')

n_d      = length(d_values_to_plot);
eig_data = zeros(n_d, length(n_wave));
valid_d  = true(n_d, 1);

for idx = 1:n_d
    d_index = find(abs(dv_range - d_values_to_plot(idx)) < 1e-9, 1);
    if isempty(d_index)
        warning('d = %g not in dv_range.', d_values_to_plot(idx))
        valid_d(idx) = false; continue
    end
    eig_data(idx,:) = all_max_eigs(d_index, :);
end

all_vals = eig_data(valid_d, :);
y_lo = min(all_vals(:)) * 1.15;
y_hi = max(all_vals(:)) * 1.15;
if y_hi < 0, y_hi =  10; end
if y_lo > 0, y_lo = -10; end

largest_d_idx = find(valid_d, 1, 'last');
curve_top     = eig_data(largest_d_idx, :);
cross_indices = find(diff(sign(curve_top)) ~= 0);

n_cross = zeros(1, length(cross_indices));
for ci = 1:length(cross_indices)
    i1 = cross_indices(ci); i2 = i1 + 1;
    n_cross(ci) = n_wave(i1) - curve_top(i1) * ...
        (n_wave(i2)-n_wave(i1)) / (curve_top(i2)-curve_top(i1));
end

all_bounds  = [n_wave(1), n_cross, n_wave(end)];
shade_color = [0.92 0.92 0.92];

for bi = 1:length(all_bounds)-1
    mid_val = interp1(n_wave, curve_top, 0.5*(all_bounds(bi)+all_bounds(bi+1)), 'linear');
    if mid_val > 0
        patch(ax, [all_bounds(bi) all_bounds(bi+1) all_bounds(bi+1) all_bounds(bi)], ...
            [y_lo y_lo y_hi y_hi], shade_color, ...
            'EdgeColor','none', 'FaceAlpha',1.0, 'HandleVisibility','off')
    end
end

plot(ax, [n_wave(1) n_wave(end)], [0 0], '--k', 'LineWidth',1.2, 'HandleVisibility','off')

line_colors = {[0.12 0.47 0.71], [1.0 0 0], [0 0.6 0]};
for idx = 1:n_d
    if ~valid_d(idx), continue; end
    plot(ax, n_wave, eig_data(idx,:), ...
        'Color',line_colors{idx}, 'LineStyle','-', 'LineWidth',2.2, ...
        'DisplayName', sprintf('$d = %g$', d_values_to_plot(idx)))
end

xlim(ax, [0 n_wave(end)]); ylim(ax, [y_lo y_hi])
xlabel(ax, '$n$',                       'Interpreter','latex', 'FontSize',16)
ylabel(ax, '$\mathrm{Re}\,\lambda(k)$', 'Interpreter','latex', 'FontSize',16)
set(ax, 'FontSize',14, 'FontName','Times New Roman', ...
    'TickLabelInterpreter','latex', 'LineWidth',1.0, 'Layer','top', ...
    'Box','on', 'XMinorTick','on', 'YMinorTick','on', 'GridAlpha',0.25)
grid(ax, 'on')
leg = legend(ax, 'Location','best', 'Interpreter','latex', 'FontSize',13);
leg.Box = 'on';

d_tag = strrep(sprintf('d_%g_%g_%g', d_values_to_plot(1), d_values_to_plot(2), d_values_to_plot(3)), '.','p');
exportgraphics(fig2, fullfile(save_folder, ...
    sprintf('%s_%s_dispersion_relation.jpg', param_tag, d_tag)), 'Resolution',400)

%% Figure 3/4 — Limit cycle plots (only if applicable)
if exist('x_limit','var') && ~isempty(x_limit)

    fig4 = figure(4);
    set(fig4, 'Color','w', 'Units','centimeters', 'Position',[2 2 18 12])
    plot(t_limit, x_limit(:,1), 'LineWidth',2.0, 'DisplayName','$R$');    hold on
    plot(t_limit, x_limit(:,2), 'LineWidth',2.0, 'DisplayName','$R_i$')
    plot(t_limit, x_limit(:,3), 'LineWidth',2.0, 'DisplayName','$\rho$')
    plot(t_limit, x_limit(:,4), 'LineWidth',2.0, 'DisplayName','$\rho_i$')
    xlabel('Time', 'Interpreter','latex', 'FontSize',14)
    ylabel('Solution', 'Interpreter','latex', 'FontSize',14)
    legend('Interpreter','latex', 'FontSize',13, 'Location','best')
    grid on; box on
    set(gca, 'FontSize',14, 'FontName','Times New Roman', ...
        'TickLabelInterpreter','latex', 'LineWidth',1.0)
    exportgraphics(fig4, fullfile(save_folder, ...
        sprintf('%s_d_NA_limit_cycle_time_series.jpg', param_tag)), 'Resolution',400)

    fig5 = figure(5);
    set(fig5, 'Color','w', 'Units','centimeters', 'Position',[2 2 14 12])
    plot(x_limit(:,1), x_limit(:,3), 'k-', 'LineWidth',2.0); hold on
    plot(x_limit(1,1), x_limit(1,3), 'ko', 'MarkerFaceColor','k', 'MarkerSize',7)
    xlabel('$R$',    'Interpreter','latex', 'FontSize',14)
    ylabel('$\rho$', 'Interpreter','latex', 'FontSize',14)
    grid on; box on
    set(gca, 'FontSize',14, 'FontName','Times New Roman', ...
        'TickLabelInterpreter','latex', 'LineWidth',1.0)
    exportgraphics(fig5, fullfile(save_folder, ...
        sprintf('%s_d_NA_limit_cycle_phase_plane_R_rho.jpg', param_tag)), 'Resolution',400)
end

%% Combined summary figure
has_lc      = exist('x_limit','var') && ~isempty(x_limit);
panelSz     = 420;
shade_comb  = [0.92 0.92 0.92];

if has_lc
    fig_comb = figure('Color','w', 'Units','pixels', 'Position',[50 50 2*panelSz 2*panelSz]);

    axA = subplot(2,2,1);
    plot(t_limit, x_limit(:,1), 'LineWidth',2, 'DisplayName','$R$');     hold on
    plot(t_limit, x_limit(:,2), 'LineWidth',2, 'DisplayName','$R_i$')
    plot(t_limit, x_limit(:,3), 'LineWidth',2, 'DisplayName','$\rho$')
    plot(t_limit, x_limit(:,4), 'LineWidth',2, 'DisplayName','$\rho_i$')
    xlabel('Time', 'Interpreter','latex', 'FontSize',13)
    ylabel('Solution', 'Interpreter','latex', 'FontSize',13)
    legend('Interpreter','latex', 'FontSize',11, 'Location','best')
    grid on; box on; pbaspect([1 1 1])
    set(axA, 'FontSize',12, 'FontName','Times New Roman', ...
        'TickLabelInterpreter','latex', 'LineWidth',0.8)

    axB = subplot(2,2,2);
    plot(x_limit(:,1), x_limit(:,3), 'k-', 'LineWidth',2); hold on
    plot(x_limit(1,1), x_limit(1,3), 'ko', 'MarkerFaceColor','k', 'MarkerSize',7)
    xlabel('$R$',    'Interpreter','latex', 'FontSize',13)
    ylabel('$\rho$', 'Interpreter','latex', 'FontSize',13)
    grid on; box on; pbaspect([1 1 1])
    set(axB, 'FontSize',12, 'FontName','Times New Roman', ...
        'TickLabelInterpreter','latex', 'LineWidth',0.8)

    axC = subplot(2,2,3);
    plot(dv_range, max_eig_vs_d, 'b-', 'LineWidth',2); hold on
    yline(0, '--k', 'LineWidth',1.2)
    if ~isempty(cross_d)
        xline(axC, d_crit, '--r', 'LineWidth',1.5)
        text(axC, d_crit+0.3, max(max_eig_vs_d)*0.9, ...
            sprintf('$\\mathbf{d_c = %.2f}$', d_crit), ...
            'Interpreter','latex', 'FontSize',14, 'Color',[0.85 0.10 0.10], ...
            'FontName','Times New Roman', 'HorizontalAlignment','left')
    end
    xlabel('Diffusion ratio, $d$',         'Interpreter','latex', 'FontSize',13)
    ylabel('$\max\,\mathrm{Re}(\lambda)$', 'Interpreter','latex', 'FontSize',13)
    grid on; box on; pbaspect([1 1 1])
    set(axC, 'FontSize',12, 'FontName','Times New Roman', ...
        'TickLabelInterpreter','latex', 'LineWidth',0.8, 'GridAlpha',0.3)

    axD = subplot(2,2,4);
    hold(axD, 'on')
    for bi = 1:length(all_bounds)-1
        mid_val = interp1(n_wave, curve_top, 0.5*(all_bounds(bi)+all_bounds(bi+1)), 'linear');
        if mid_val > 0
            patch(axD, [all_bounds(bi) all_bounds(bi+1) all_bounds(bi+1) all_bounds(bi)], ...
                [y_lo y_lo y_hi y_hi], shade_comb, ...
                'EdgeColor','none', 'FaceAlpha',1.0, 'HandleVisibility','off')
        end
    end
    plot(axD, [n_wave(1) n_wave(end)],[0 0], '--k', 'LineWidth',1.2, 'HandleVisibility','off')
    for idx = 1:n_d
        if ~valid_d(idx), continue; end
        plot(axD, n_wave, eig_data(idx,:), ...
            'Color',line_colors{idx}, 'LineStyle','-', 'LineWidth',2.2, ...
            'DisplayName', sprintf('$d = %g$', d_values_to_plot(idx)))
    end
    xlim(axD,[0 n_wave(end)]); ylim(axD,[y_lo y_hi])
    xlabel(axD,'$n$',                       'Interpreter','latex', 'FontSize',13)
    ylabel(axD,'$\mathrm{Re}\,\lambda(k)$', 'Interpreter','latex', 'FontSize',13)
    legend(axD, 'Location','best', 'Interpreter','latex', 'FontSize',11, 'Box','on')
    grid(axD,'on'); pbaspect([1 1 1])
    set(axD, 'FontSize',12, 'FontName','Times New Roman', ...
        'TickLabelInterpreter','latex', 'LineWidth',0.8, 'GridAlpha',0.25, ...
        'Layer','top', 'Box','on', 'XMinorTick','on', 'YMinorTick','on')

    drawnow
    addPanelLabel(axA,'A',fig_comb); addPanelLabel(axB,'B',fig_comb)
    addPanelLabel(axC,'C',fig_comb); addPanelLabel(axD,'D',fig_comb)

else
    fig_comb = figure('Color','w', 'Units','pixels', 'Position',[50 50 2*panelSz panelSz]);

    axA = subplot(1,2,1);
    plot(dv_range, max_eig_vs_d, 'b-', 'LineWidth',2); hold on
    yline(0, '--k', 'LineWidth',1.2)
    if ~isempty(cross_d)
        xline(axA, d_crit, '--r', 'LineWidth',1.5)
        text(axA, d_crit+0.3, max(max_eig_vs_d)*0.9, ...
            sprintf('$\\mathbf{d_c = %.2f}$', d_crit), ...
            'Interpreter','latex', 'FontSize',14, 'Color',[1 0 0], ...
            'FontName','Times New Roman', 'HorizontalAlignment','left')
    end
    xlabel('Diffusion ratio, $d$',         'Interpreter','latex', 'FontSize',13)
    ylabel('$\max\,\mathrm{Re}(\lambda)$', 'Interpreter','latex', 'FontSize',13)
    grid on; box on; pbaspect([1 1 1])
    set(axA, 'FontSize',12, 'FontName','Times New Roman', ...
        'TickLabelInterpreter','latex', 'LineWidth',0.8, 'GridAlpha',0.3)

    axB = subplot(1,2,2);
    hold(axB, 'on')
    for bi = 1:length(all_bounds)-1
        mid_val = interp1(n_wave, curve_top, 0.5*(all_bounds(bi)+all_bounds(bi+1)), 'linear');
        if mid_val > 0
            patch(axB, [all_bounds(bi) all_bounds(bi+1) all_bounds(bi+1) all_bounds(bi)], ...
                [y_lo y_lo y_hi y_hi], shade_comb, ...
                'EdgeColor','none', 'FaceAlpha',1.0, 'HandleVisibility','off')
        end
    end
    plot(axB, [n_wave(1) n_wave(end)],[0 0], '--k', 'LineWidth',1.2, 'HandleVisibility','off')
    for idx = 1:n_d
        if ~valid_d(idx), continue; end
        plot(axB, n_wave, eig_data(idx,:), ...
            'Color',line_colors{idx}, 'LineStyle','-', 'LineWidth',2.2, ...
            'DisplayName', sprintf('$d = %g$', d_values_to_plot(idx)))
    end
    xlim(axB,[0 n_wave(end)]); ylim(axB,[y_lo y_hi])
    xlabel(axB,'$n$',                       'Interpreter','latex', 'FontSize',13)
    ylabel(axB,'$\mathrm{Re}\,\lambda(k)$', 'Interpreter','latex', 'FontSize',13)
    legend(axB, 'Location','best', 'Interpreter','latex', 'FontSize',11, 'Box','on')
    grid(axB,'on'); pbaspect([1 1 1])
    set(axB, 'FontSize',12, 'FontName','Times New Roman', ...
        'TickLabelInterpreter','latex', 'LineWidth',0.8, 'GridAlpha',0.25, ...
        'Layer','top', 'Box','on', 'XMinorTick','on', 'YMinorTick','on')

    drawnow
    addPanelLabel(axA,'A',fig_comb); addPanelLabel(axB,'B',fig_comb)
end

exportgraphics(fig_comb, fullfile(save_folder, ...
    sprintf('%s_combined_summary.jpg', param_tag)), 'Resolution',400)


%% -----------------------------------------------------------------------
%  Local functions
%% -----------------------------------------------------------------------

function res = rho_kinetics_reduced(x, p)
    R = x(1);  rho = x(2);
    Ri = 1-R;  rhoi = 1-rho;
    n = p.n;
    f = p.a1*Ri^n/(p.k1^n+Ri^n) + p.a2*Ri^n*R/(p.k2^n+Ri^n) ...
        - p.a3*rho*R^n/(p.k3^n+R^n) - R;
    g = p.a4*rhoi^n/(p.k4^n+rhoi^n) + p.a5*R*rhoi^n/(p.k5^n+rhoi^n) - p.a6*rho;
    res = [f; g];
end

function dxdt = Nondim_full(~, x, p, gamma)
    R=x(1); Ri=x(2); rho=x(3); rhoi=x(4);
    n = p.n;
    f = p.a1*Ri^n/(p.k1^n+Ri^n) + p.a2*Ri^n*R/(p.k2^n+Ri^n) ...
        - p.a3*rho*R^n/(p.k3^n+R^n) - R;
    g = p.a4*rhoi^n/(p.k4^n+rhoi^n) + p.a5*R*rhoi^n/(p.k5^n+rhoi^n) - p.a6*rho;
    dxdt = gamma*[f; -f; g; -g];
end

function J = get_4x4_jacobian_full(R, Ri, rho, rhoi, p, gamma)
    n = p.n;
    f_R    =  p.a2*Ri^n/(p.k2^n+Ri^n) - p.a3*rho*n*p.k3^n*R^(n-1)/(p.k3^n+R^n)^2 - 1;
    f_Ri   =  p.a1*n*p.k1^n*Ri^(n-1)/(p.k1^n+Ri^n)^2 + p.a2*R*n*p.k2^n*Ri^(n-1)/(p.k2^n+Ri^n)^2;
    f_rho  = -p.a3*R^n/(p.k3^n+R^n);
    f_rhoi =  0;
    g_R    =  p.a5*rhoi^n/(p.k5^n+rhoi^n);
    g_Ri   =  0;
    g_rho  = -p.a6;
    g_rhoi =  p.a4*n*p.k4^n*rhoi^(n-1)/(p.k4^n+rhoi^n)^2 + p.a5*R*n*p.k5^n*rhoi^(n-1)/(p.k5^n+rhoi^n)^2;
    J = gamma * [ f_R,  f_Ri,  f_rho,  f_rhoi;
                 -f_R, -f_Ri, -f_rho, -f_rhoi;
                  g_R,  g_Ri,  g_rho,  g_rhoi;
                 -g_R, -g_Ri, -g_rho, -g_rhoi];
end

function Jred = get_reduced_jacobian_numeric(x, p, gamma)
    h = 1e-6; m = length(x); Jred = zeros(m,m);
    for j = 1:m
        e = zeros(m,1); e(j) = 1;
        Jred(:,j) = (gamma*rho_kinetics_reduced(x+h*e,p) - ...
                     gamma*rho_kinetics_reduced(x-h*e,p)) / (2*h);
    end
end

function [x_limit, t_limit, T_final, deviation] = compute_limit_cycle_4var(t, x)
    if any(x(:) < 0)
        x_limit=[]; t_limit=[]; T_final=NaN; deviation=NaN; return
    end
    [~,pk1] = findpeaks(x(:,1), 'MinPeakDistance',5); pt1 = t(pk1);
    [~,pk2] = findpeaks(x(:,3), 'MinPeakDistance',5); pt2 = t(pk2);
    if length(pt1) < 3 || length(pt2) < 3
        x_limit=[]; t_limit=[]; T_final=NaN; deviation=NaN; return
    end
    T1 = mean(diff(pt1(2:end)));
    T2 = mean(diff(pt2(2:end)));
    is1=find(t>=pt1(end-1),1); ie1=find(t<=pt1(end),1,'last');
    xl1=x(is1:ie1,:); tl1=t(is1:ie1)-t(is1);
    is2=find(t>=pt2(end-1),1); ie2=find(t<=pt2(end),1,'last');
    xl2=x(is2:ie2,:); tl2=t(is2:ie2)-t(is2);
    d1=norm(xl1(1,:)-xl1(end,:)); d2=norm(xl2(1,:)-xl2(end,:));
    if d1 < d2
        x_limit=xl1; t_limit=tl1; T_final=T1; deviation=d1;
    else
        x_limit=xl2; t_limit=tl2; T_final=T2; deviation=d2;
    end
end

function Javg = average_jacobian_on_cycle(t_limit, x_limit, p, gamma)
    N = length(t_limit);
    Jstore = zeros(4,4,N);
    for q = 1:N
        Jstore(:,:,q) = get_4x4_jacobian_full( ...
            x_limit(q,1), x_limit(q,2), x_limit(q,3), x_limit(q,4), p, gamma);
    end
    T = t_limit(end) - t_limit(1);
    Javg = zeros(4,4);
    for a = 1:4
        for b = 1:4
            Javg(a,b) = trapz(t_limit, squeeze(Jstore(a,b,:))) / T;
        end
    end
end

function addPanelLabel(ax, lbl, fig)
    pos = get(ax, 'Position');
    annotation(fig, 'textbox', [pos(1), pos(2)+pos(4), 0.05, 0.05], ...
        'String',              ['\bf' lbl], ...
        'FontName',            'Times New Roman', ...
        'FontSize',            14, ...
        'FontWeight',          'bold', ...
        'EdgeColor',           'none', ...
        'BackgroundColor',     'none', ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment',   'bottom', ...
        'Interpreter',         'tex');
end
