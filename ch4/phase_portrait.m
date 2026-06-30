%% Multi-scenario panel: S1, S2, S3
% Plots phase portraits and time series for the Rac/Rho ODE system
% across the three dynamical regimes (bistable, oscillatory, coexistence)

clear; clc; close all;

scenario = 'S3'; % switch to 'S1' or 'S2' to see the other regimes
gamma = 1;
col_R   = [0.85 0.15 0.10];
col_rho = [0.15 0.40 0.85];
t_span  = [0 15];

folderName = fullfile(pwd, 'Nullcline_Analysis_ODE');

if strcmp(scenario, 'S1')
    a1 = 0.45; k1 = 1; a2 = 6; k2 = 0.02; n = 1;
    a4 = 1.6; k4 = 1.5; k5 = 0.75; a3 = 5.5; k3 = 0.3; a6 = 1;
    a5_values = [1, 5, 15];
elseif strcmp(scenario, 'S2')
    a1 = 0.25; k1 = 2.5; a2 = 5.5; k2 = 0.075; n = 1;
    a4 = 0.1; k4 = 1.5; k5 = 0.55; a3 = 6; k3 = 0.3; a6 = 1;
    a5_values = [0.5, 2.3, 5];
elseif strcmp(scenario, 'S3')
    a1 = 0.5; k1 = 3; a2 = 7; k2 = 0.075; n = 1;
    a4 = 0.2; k4 = 2; k5 = 0.5; a3 = 4.7; k3 = 0.25; a6 = 0.6;
    a5_values = [1, 4.3, 12];
else
    error('Invalid scenario.');
end

fig = figure('Color','w','Position',[100 50 900 1200]);

for i = 1:length(a5_values)
    a5 = a5_values(i);

    g = @(R, rho) gamma * ( ...
        (a4 .* (1-rho).^n) ./ (k4.^n + (1-rho).^n) + ...
        (a5 .* R .* (1-rho).^n) ./ (k5.^n + (1-rho).^n) - ...
        a6 .* rho );
    f = @(R, rho) gamma * ( ...
        (a1 .* (1-R).^n) ./ (k1.^n + (1-R).^n) + ...
        (a2 .* (1-R).^n .* R) ./ (k2.^n + (1-R).^n) - ...
        (a3 .* rho .* R.^n) ./ (k3.^n + R.^n) - R );
    ode_sys = @(t, y) [f(y(1), y(2)); g(y(1), y(2))];

    % use two initial conditions for the middle a5 value in S1/S3 to show bistability
    if ismember(scenario, {'S3', 'S1'}) && i == 2
        IC_list = [0.1 0.7; 0.7 0.6];
    else
        IC_list = [0.1 0.1];
    end

    % phase portrait
    subplot(length(a5_values), 2, 2*i-1); hold on;
    [R_grid, rho_grid] = meshgrid(linspace(0,1,25));
    dR = f(R_grid, rho_grid); drho = g(R_grid, rho_grid);
    vel = sqrt(dR.^2 + drho.^2); vel(vel==0) = 1;
    quiver(R_grid, rho_grid, dR./vel, drho./vel, 0.5, 'Color',[0.7 0.7 0.7]);

    [R_fg, rho_fg] = meshgrid(linspace(0,1,300));
    contour(R_fg, rho_fg, f(R_fg, rho_fg), [0 0], 'Color', col_R,   'LineWidth', 2.5);
    contour(R_fg, rho_fg, g(R_fg, rho_fg), [0 0], 'Color', col_rho, 'LineWidth', 2.5);

    for j = 1:size(IC_list,1)
        [~, Y_sol_j] = ode15s(ode_sys, t_span, IC_list(j,:));
        plot(Y_sol_j(:,1), Y_sol_j(:,2), 'k', 'LineWidth', 2);
        plot(IC_list(j,1), IC_list(j,2), 'o', 'MarkerFaceColor','g', 'MarkerSize', 7);
    end
    xlabel('$R$', 'Interpreter', 'latex');
    ylabel('$\rho$', 'Interpreter', 'latex');
    axis([0 1 0 1]); box on;

    % time series
    subplot(length(a5_values), 2, 2*i); hold on;
    legend_entries = {};

    for j = 1:size(IC_list,1)
        [t_sol, Y_sol] = ode15s(ode_sys, t_span, IC_list(j,:));

        if j == 1
            lt = '--';
        else
            lt = '-';
        end

        plot(t_sol, Y_sol(:,1), 'Color', col_R,   'LineWidth', 2, 'LineStyle', lt);
        plot(t_sol, Y_sol(:,2), 'Color', col_rho, 'LineWidth', 2, 'LineStyle', lt);

        legend_entries{end+1} = sprintf('$R(t)$ IC%d', j);
        legend_entries{end+1} = sprintf('$\\rho(t)$ IC%d', j);
    end

    xlabel('Time', 'Interpreter', 'latex');
    ylabel('Concentration', 'Interpreter', 'latex');
    legend(legend_entries, 'Interpreter', 'latex', 'Location', 'best');
    xlim(t_span);
    box on;
end

sgtitle(['Scenario ' scenario ': $a_5 = ' sprintf('%g, ', a5_values(1:end-1)) ...
         sprintf('%g$', a5_values(end))], 'Interpreter', 'latex', 'FontSize', 14);

if ~exist(folderName, 'dir'), mkdir(folderName); end
timestamp = datestr(now,'yyyymmdd_HHMMSS');
savePath = fullfile(folderName, [scenario '_Panel_' timestamp '.png']);

exportgraphics(fig, savePath, 'Resolution', 400);
fprintf('Saved figure to: %s\n', savePath);
