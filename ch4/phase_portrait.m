%% MULTI-SCENARIO PANEL: S1, S2, S3
clear; clc; close all;

% --- Configuration ---
scenario = 'S3'; % Change to 'S1' or 'S2' as needed
pp.gamma = 1;
col_R   = [0.85 0.15 0.10];
col_rho = [0.15 0.40 0.85];
t_span  = [0 15]; % Increased range as requested (1-15)

% Define target directory
folderName = fullfile(pwd, 'Nullcline_Analysis_ODE');

switch scenario
    case 'S1'
        base.a1 = 0.45; base.k1 = 1; base.a2 = 6;
        base.k2 = 0.02; base.n = 1; base.a4 = 1.6;
        base.k4 = 1.5; base.k5 = 0.75; base.a3 = 5.5;
        base.k3 = 0.3; base.a6 = 1;
        a5_values = [1, 5, 15];
    case 'S2'
        base.a1 = 0.25; base.k1 = 2.5; base.a2 = 5.5;
        base.k2 = 0.075; base.n = 1; base.a4 = 0.1;
        base.k4 = 1.5; base.k5 = 0.55; base.a3 = 6;
        base.k3 = 0.3; base.a6 = 1;
        a5_values = [0.5, 2.3, 5];
    case 'S3'
        base.a1 = 0.5; base.k1 = 3; base.a2 = 7;
        base.k2 = 0.075; base.n = 1; base.a4 = 0.2;
        base.k4 = 2; base.k5 = 0.5; base.a3 = 4.7;
        base.k3 = 0.25; base.a6 = 0.6;
        a5_values = [1, 4.3, 12];
    otherwise
        error('Invalid scenario.');
end

fig = figure('Color','w','Position',[100 50 900 1200]);

for i = 1:length(a5_values)
    pp_current = base;
    pp_current.a5 = a5_values(i);
    pp_current.gamma = pp.gamma;
   
    % System dynamics
    g = @(R, rho) pp_current.gamma * ( ...
        (pp_current.a4 .* (1-rho).^pp_current.n) ./ (pp_current.k4.^pp_current.n + (1-rho).^pp_current.n) + ...
        (pp_current.a5 .* R .* (1-rho).^pp_current.n) ./ (pp_current.k5.^pp_current.n + (1-rho).^pp_current.n) - ...
        pp_current.a6 .* rho );
    f = @(R, rho) pp_current.gamma * ( ...
        (pp_current.a1 .* (1-R).^pp_current.n) ./ (pp_current.k1.^pp_current.n + (1-R).^pp_current.n) + ...
        (pp_current.a2 .* (1-R).^pp_current.n .* R) ./ (pp_current.k2.^pp_current.n + (1-R).^pp_current.n) - ...
        (pp_current.a3 .* rho .* R.^pp_current.n) ./ (pp_current.k3.^pp_current.n + R.^pp_current.n) - R );
    ode_sys = @(t, y) [f(y(1), y(2)); g(y(1), y(2))];
   
    % Initial Condition selection
    if ismember(scenario, {'S3', 'S1'}) && i == 2
        IC_list = [0.1 0.7; 0.7 0.6]; % Lower IC (1) and Upper IC (2)
    else
        IC_list = [0.1 0.1];
    end
    
    % --- Phase Portrait ---
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
    
    % --- Time Series ---
    subplot(length(a5_values), 2, 2*i); hold on;
    legend_entries = {};
    
    for j = 1:size(IC_list,1)
        [t_sol, Y_sol] = ode15s(ode_sys, t_span, IC_list(j,:));
        
        % Logic: Dashed (--) for lower IC1, Solid (-) for upper IC2
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
    xlim(t_span); % Set x-axis limit to 15
    box on;
end

sgtitle(['Scenario ' scenario ': $a_5 = ' sprintf('%g, ', a5_values(1:end-1)) ...
         sprintf('%g$', a5_values(end))], 'Interpreter', 'latex', 'FontSize', 14);

% --- Saving Process ---
if ~exist(folderName, 'dir'), mkdir(folderName); end
timestamp = datestr(now,'yyyymmdd_HHMMSS');
savePath = fullfile(folderName, [scenario '_Panel_' timestamp '.png']);

exportgraphics(fig, savePath, 'Resolution', 400);
fprintf('Saved figure to: %s\n', savePath);