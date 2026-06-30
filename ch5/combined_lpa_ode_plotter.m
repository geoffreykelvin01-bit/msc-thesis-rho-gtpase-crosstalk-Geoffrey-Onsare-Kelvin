%% plot matcont files for LPA and ODE - Full Professional Suite
clc; clear; close all
% =========================================================================
% 0. CONFIGURATION & FILE PATHS
% =========================================================================
set_num = 'set_b';
file_suffix = ['_', set_num]; 

baseDir = fullfile(pwd, 'OSCILLATORY_trial');

fig_folder = fullfile(baseDir, ['figures_for_', set_num]);
if ~exist(fig_folder, 'dir'), mkdir(fig_folder); end

% Axis Limits
lim.a5     = [0, 13]; 
lim.a52    = [0, 16]; 
lim.a2     = [2, 14];
lim.y_ode  = [-0.1, 1.5];     
lim.y_lpa  = [-0.1, 1.5];   

% Standard Labels
titles_ode = {'R', '\rho'}; 
titles_lpa = {'R_l', '\rho_l'}; 
file_names_vars = {'Rac', 'Rho'}; 

pathLPA = fullfile(baseDir, ['LPA', file_suffix, '\']);
pathODE = fullfile(baseDir, ['ODE', file_suffix, '\']);

% =========================================================================
% ROW INDEX REFERENCE
%   ODE 1-param:  row1=R,  row2=rho,  row3=a5
%   ODE 2-param:  row1=R,  row2=rho,  row3=a2,  row4=a5
%   LPA 1-param:  row1=Rl, row2=rho_l, row3=Rg, row4=rho_g, row5=a5
%   LPA 2-param:  row1=Rl, row2=rho_l, row3=Rg, row4=rho_g, row5=a2, row6=a5
% =========================================================================
ODE1_A5_ROW  = 3;
ODE2_A2_ROW  = 3;  ODE2_A5_ROW  = 4;
LPA1_A5_ROW  = 5;
LPA2_A2_ROW  = 5;  LPA2_A5_ROW  = 6;

LPA1_RG_ROW   = 3;
LPA1_RHOG_ROW = 4;

% =========================================================================
% 1. DATA LOADING
% =========================================================================
ode_ep       = combineFiles(pathODE, 'EP_EP');
ode_hh       = combineHH(pathODE, 'H_H', ODE2_A5_ROW, ODE2_A2_ROW);
ode_lp_segs  = loadIndependently(pathODE, 'LP_LP');
ode_lp_lp    = combineAllFiles(pathODE, 'LP_LP');
ode_lc_cells = loadIndependently(pathODE, 'H_LC');

lpa_ep       = combineFiles(pathLPA, 'EP_EP');
lpa_bp       = combineFiles(pathLPA, 'BP_EP');
lpa_hh       = combineHH(pathLPA, 'H_H', LPA2_A5_ROW, LPA2_A2_ROW);
lpa_lp_segs  = loadIndependently(pathLPA, 'LP_LP');
lpa_lp_lp    = combineAllFiles(pathLPA, 'LP_LP');
lpa_lc_cells = loadIndependently(pathLPA, 'H_LC');

% Diagnostic printout
fprintf('ODE H_H points: %d\n', size(ode_hh.x, 2));
fprintf('LPA H_H points: %d\n', size(lpa_hh.x, 2));

% =========================================================================
% 2. FIGURES 1 & 2: ODE 1-parameter
% =========================================================================
for p = 1:2
    fig = figure('Color','w','Units','pixels','Position',[100 100 800 600]);
    hold on; box on; set(gca,'Color','w','FontSize',12);
    plotBranch(ode_ep, p, ODE1_A5_ROW, 'b', 'r');
    plotLC_Envelope(ode_lc_cells, p, 2, 'g');
    labelBifurcations(ode_ep, p, ODE1_A5_ROW);
    xlabel('$a_5$','Interpreter','latex','FontSize',16);
    ylabel(titles_ode{p},'FontSize',16);
    xlim(lim.a5); ylim(lim.y_ode); grid on;
    exportgraphics(fig, fullfile(fig_folder, sprintf('ODE_%s%s.jpg',file_names_vars{p},file_suffix)), 'Resolution',300);
end

% =========================================================================
% 3. FIGURES 3 & 4: LPA 1-parameter
% =========================================================================
for p = 1:2
    fig = figure('Color','w','Units','pixels','Position',[100 100 800 600]);
    hold on; box on; set(gca,'Color','w','FontSize',12);
    plotBranch(lpa_ep, p, LPA1_A5_ROW, 'b', 'r');
    plotBranch(lpa_bp, p, LPA1_A5_ROW, 'c', 'k');
    plotLC_Envelope(lpa_lc_cells, p, 4, 'g');
    labelBifurcations(lpa_ep, p, LPA1_A5_ROW);
    xlabel('$a_5$','Interpreter','latex','FontSize',16);
    ylabel(titles_lpa{p},'FontSize',16);
    xlim(lim.a5); ylim(lim.y_lpa); grid on;
    exportgraphics(fig, fullfile(fig_folder, sprintf('LPA_%s%s.jpg',file_names_vars{p},file_suffix)), 'Resolution',300);
end

% =========================================================================
% 4. PANEL FIGURE (4-panel 2x2): (a) R_l  (b) rho_l  (c) R  (d) rho
% =========================================================================
panel_labels = {'A','B','C','D'};
figPanel = figure('Color','w','Units','inches',...
    'Position',[0 0 18 14],'PaperUnits','inches','PaperSize',[18 14],'PaperPosition',[0 0 18 14]);
tl = tiledlayout(figPanel,2,2,'TileSpacing','compact','Padding','loose');

ax1 = nexttile(tl); hold on; box on; set(ax1,'Color','w','FontSize',12);
plotBranch(lpa_ep,1,LPA1_A5_ROW,'b','r');
plotBranch(lpa_bp,1,LPA1_A5_ROW,'c','k');
plotLC_Envelope(lpa_lc_cells,1,4,'g');
labelBifurcations(lpa_ep,1,LPA1_A5_ROW);
xlabel('$a_5$','Interpreter','latex','FontSize',16);
ylabel('$R_l$','Interpreter','latex','FontSize',16);
xlim(lim.a5); ylim(lim.y_lpa); grid on;
addPanelLabel(figPanel, ax1, panel_labels{1});

ax2 = nexttile(tl); hold on; box on; set(ax2,'Color','w','FontSize',12);
plotBranch(lpa_ep,2,LPA1_A5_ROW,'b','r');
plotBranch(lpa_bp,2,LPA1_A5_ROW,'c','k');
plotLC_Envelope(lpa_lc_cells,2,4,'g');
labelBifurcations(lpa_ep,2,LPA1_A5_ROW);
xlabel('$a_5$','Interpreter','latex','FontSize',16);
ylabel('$\rho_l$','Interpreter','latex','FontSize',16);
xlim(lim.a5); ylim(lim.y_lpa); grid on;
addPanelLabel(figPanel, ax2, panel_labels{2});

ax3 = nexttile(tl); hold on; box on; set(ax3,'Color','w','FontSize',12);
plotBranch(ode_ep,1,ODE1_A5_ROW,'b','r');
plotLC_Envelope(ode_lc_cells,1,2,'g');
labelBifurcations(ode_ep,1,ODE1_A5_ROW);
xlabel('$a_5$','Interpreter','latex','FontSize',16);
ylabel('$R$','Interpreter','latex','FontSize',16);
xlim(lim.a5); ylim(lim.y_ode); grid on;
addPanelLabel(figPanel, ax3, panel_labels{3});

ax4 = nexttile(tl); hold on; box on; set(ax4,'Color','w','FontSize',12);
plotBranch(ode_ep,2,ODE1_A5_ROW,'b','r');
plotLC_Envelope(ode_lc_cells,2,2,'g');
labelBifurcations(ode_ep,2,ODE1_A5_ROW);
xlabel('$a_5$','Interpreter','latex','FontSize',16);
ylabel('$\rho$','Interpreter','latex','FontSize',16);
xlim(lim.a5); ylim(lim.y_ode); grid on;
addPanelLabel(figPanel, ax4, panel_labels{4});

exportgraphics(figPanel, fullfile(fig_folder,['PanelFigure_ABCD',file_suffix,'.jpg']),'Resolution',300);

% =========================================================================
% 4b. PANEL FIGURE (2-panel): (a) R_l  (b) R
% =========================================================================
figPanel2 = figure('Color','w','Units','inches',...
    'Position',[0 0 18 7],'PaperUnits','inches','PaperSize',[18 7],'PaperPosition',[0 0 18 7]);
tl2 = tiledlayout(figPanel2,1,2,'TileSpacing','compact','Padding','loose');

ax2a = nexttile(tl2); hold on; box on; set(ax2a,'Color','w','FontSize',12);
plotBranch(lpa_ep,1,LPA1_A5_ROW,'b','r');
plotBranch(lpa_bp,1,LPA1_A5_ROW,'c','k');
plotLC_Envelope(lpa_lc_cells,1,4,'g');
labelBifurcations(lpa_ep,1,LPA1_A5_ROW);
xlabel('$a_5$','Interpreter','latex','FontSize',16);
ylabel('$R_l$','Interpreter','latex','FontSize',16);
xlim(lim.a5); ylim(lim.y_lpa); grid on;
addPanelLabel(figPanel2, ax2a, 'A');

ax2b = nexttile(tl2); hold on; box on; set(ax2b,'Color','w','FontSize',12);
plotBranch(ode_ep,1,ODE1_A5_ROW,'b','r');
plotLC_Envelope(ode_lc_cells,1,2,'g');
labelBifurcations(ode_ep,1,ODE1_A5_ROW);
xlabel('$a_5$','Interpreter','latex','FontSize',16);
ylabel('$R$','Interpreter','latex','FontSize',16);
xlim(lim.a5); ylim(lim.y_ode); grid on;
addPanelLabel(figPanel2, ax2b, 'B');

exportgraphics(figPanel2, fullfile(fig_folder,['PanelFigure_AB',file_suffix,'.jpg']),'Resolution',300);

% =========================================================================
% 5. FIGURE 5: TWO-PARAMETER OVERLAY (HH only)
% =========================================================================
fig5 = figure('Color','w','Units','pixels','Position',[100 100 800 600]);
hold on; box on; set(gca,'Color','w','FontSize',12);
if ~isempty(ode_hh.x)
    x = ode_hh.x(ODE2_A5_ROW,:);
    y = ode_hh.x(ODE2_A2_ROW,:);
    v = isfinite(x) & isfinite(y);
    x = x(v); y = y(v);
    fill([x, fliplr(x), x(1)], [y, fliplr(y), y(1)], ...
        blendWithWhite([0 0 1],0.5), 'EdgeColor','none','FaceAlpha',1,'HandleVisibility','off');
    plot(x, y, 'b-', 'LineWidth', 2, 'HandleVisibility','off');
end
if ~isempty(lpa_hh.x)
    x = lpa_hh.x(LPA2_A5_ROW,:);
    y = lpa_hh.x(LPA2_A2_ROW,:);
    v = isfinite(x) & isfinite(y);
    x = x(v); y = y(v);
    fill([x, fliplr(x), x(1)], [y, fliplr(y), y(1)], ...
        blendWithWhite([1 0 0],0.5), 'EdgeColor','none','FaceAlpha',1,'HandleVisibility','off');
    plot(x, y, 'r-', 'LineWidth', 2, 'HandleVisibility','off');
end
xlabel('$a_5$','Interpreter','latex','FontSize',16);
ylabel('$a_2$','Interpreter','latex','FontSize',16);
grid on; xlim(lim.a52); ylim(lim.a2);
exportgraphics(fig5, fullfile(fig_folder,['TwoParam_Overlay',file_suffix,'.jpg']),'Resolution',300);

% =========================================================================
% 6. FIGURE 6: ODE 2-PARAM (LP + HH)
% =========================================================================
fig6 = figure('Color','w','Units','pixels','Position',[100 100 800 600]);
hold on; box on; set(gca,'Color','w','FontSize',12);
plotTwoParamODE(ode_hh, ode_lp_segs, lim, ODE2_A5_ROW, ODE2_A2_ROW);
xlabel('$a_5$','Interpreter','latex','FontSize',16);
ylabel('$a_2$','Interpreter','latex','FontSize',16);
grid on; xlim(lim.a52); ylim(lim.a2);
exportgraphics(fig6, fullfile(fig_folder,['ODE_2Param_Standalone',file_suffix,'.jpg']),'Resolution',300);

% =========================================================================
% 7. FIGURE 7: LPA 2-PARAM (LP + HH)
% =========================================================================
fig7 = figure('Color','w','Units','pixels','Position',[100 100 800 600]);
hold on; box on; set(gca,'Color','w','FontSize',12);
plotTwoParamLPA(lpa_hh, lpa_lp_segs, lim, LPA2_A5_ROW, LPA2_A2_ROW);
xlabel('$a_5$','Interpreter','latex','FontSize',16);
ylabel('$a_2$','Interpreter','latex','FontSize',16);
grid on; xlim(lim.a52); ylim(lim.a2);
exportgraphics(fig7, fullfile(fig_folder,['LPA_2Param_Standalone',file_suffix,'.jpg']),'Resolution',300);

% =========================================================================
% 8a. ODE THREE-PANEL: (a) R  (b) rho  (c) 2-param
% =========================================================================
figODE3 = figure('Color','w','Units','inches',...
    'Position',[0 0 24 8],'PaperUnits','inches','PaperSize',[24 8],'PaperPosition',[0 0 24 8]);
tl3 = tiledlayout(figODE3,1,3,'TileSpacing','compact','Padding','loose');

ax_o1 = nexttile(tl3); hold on; box on; set(ax_o1,'Color','w','FontSize',12);
plotBranch(ode_ep,1,ODE1_A5_ROW,'b','r'); plotLC_Envelope(ode_lc_cells,1,2,'g');
labelBifurcations(ode_ep,1,ODE1_A5_ROW);
xlabel('$a_5$','Interpreter','latex','FontSize',16);
ylabel('$R$','Interpreter','latex','FontSize',16);
xlim(lim.a5); ylim(lim.y_ode); grid on; pbaspect([1 1 1]);
addPanelLabel(figODE3, ax_o1, 'A');

ax_o2 = nexttile(tl3); hold on; box on; set(ax_o2,'Color','w','FontSize',12);
plotBranch(ode_ep,2,ODE1_A5_ROW,'b','r'); plotLC_Envelope(ode_lc_cells,2,2,'g');
labelBifurcations(ode_ep,2,ODE1_A5_ROW);
xlabel('$a_5$','Interpreter','latex','FontSize',16);
ylabel('$\rho$','Interpreter','latex','FontSize',16);
xlim(lim.a5); ylim(lim.y_ode); grid on; pbaspect([1 1 1]);
addPanelLabel(figODE3, ax_o2, 'B');

ax_o3 = nexttile(tl3); hold on; box on; set(ax_o3,'Color','w','FontSize',12);
plotTwoParamODE(ode_hh, ode_lp_segs, lim, ODE2_A5_ROW, ODE2_A2_ROW);
xlabel('$a_5$','Interpreter','latex','FontSize',16);
ylabel('$a_2$','Interpreter','latex','FontSize',16);
xlim(lim.a52); ylim(lim.a2); grid on; pbaspect([1 1 1]);
addPanelLabel(figODE3, ax_o3, 'C');

exportgraphics(figODE3, fullfile(fig_folder,['ODE_ThreePanel_abc',file_suffix,'.jpg']),'Resolution',300);

% =========================================================================
% 8b. LPA THREE-PANEL: (a) R  (b) rho  (c) 2-param
% =========================================================================
figLPA3 = figure('Color','w','Units','inches',...
    'Position',[0 0 24 8],'PaperUnits','inches','PaperSize',[24 8],'PaperPosition',[0 0 24 8]);
tl4 = tiledlayout(figLPA3,1,3,'TileSpacing','compact','Padding','loose');

ax_l1 = nexttile(tl4); hold on; box on; set(ax_l1,'Color','w','FontSize',12);
plotBranch(lpa_ep,1,LPA1_A5_ROW,'b','r');
plotBranch(lpa_bp,1,LPA1_A5_ROW,'c','k');
plotLC_Envelope(lpa_lc_cells,1,4,'g');
labelBifurcations(lpa_ep,1,LPA1_A5_ROW);
xlabel('$a_5$','Interpreter','latex','FontSize',16);
ylabel('$R$','Interpreter','latex','FontSize',16);
xlim(lim.a5); ylim(lim.y_lpa); grid on; pbaspect([1 1 1]);
addPanelLabel(figLPA3, ax_l1, 'A');

ax_l2 = nexttile(tl4); hold on; box on; set(ax_l2,'Color','w','FontSize',12);
plotBranch(lpa_ep,2,LPA1_A5_ROW,'b','r');
plotBranch(lpa_bp,2,LPA1_A5_ROW,'c','k');
plotLC_Envelope(lpa_lc_cells,2,4,'g');
labelBifurcations(lpa_ep,2,LPA1_A5_ROW);
xlabel('$a_5$','Interpreter','latex','FontSize',16);
ylabel('$\rho$','Interpreter','latex','FontSize',16);
xlim(lim.a5); ylim(lim.y_lpa); grid on; pbaspect([1 1 1]);
addPanelLabel(figLPA3, ax_l2, 'B');

ax_l3 = nexttile(tl4); hold on; box on; set(ax_l3,'Color','w','FontSize',12);
plotTwoParamOverlay(ode_hh,ode_lp_segs,lpa_hh,lpa_lp_segs,lim,...
    ODE2_A5_ROW,ODE2_A2_ROW,LPA2_A5_ROW,LPA2_A2_ROW);
xlabel('$a_5$','Interpreter','latex','FontSize',16);
ylabel('$a_2$','Interpreter','latex','FontSize',16);
xlim(lim.a52); ylim(lim.a2); grid on; pbaspect([1 1 1]);
addPanelLabel(figLPA3, ax_l3, 'C');

exportgraphics(figLPA3, fullfile(fig_folder,['LPA_ThreePanel_abc',file_suffix,'.jpg']),'Resolution',300);

% =========================================================================
% 9. ODE+LPA TWO-PARAMETER 3-PANEL: (a) ODE  (b) LPA  (c) Overlay
% =========================================================================
figTP2 = figure('Color','w','Units','inches',...
    'Position',[0 0 24 8],'PaperUnits','inches','PaperSize',[24 8],'PaperPosition',[0 0 24 8]);
tl6 = tiledlayout(figTP2,1,3,'TileSpacing','compact','Padding','loose');

ax_t1 = nexttile(tl6); hold on; box on; set(ax_t1,'Color','w','FontSize',12);
plotTwoParamODE(ode_hh, ode_lp_segs, lim, ODE2_A5_ROW, ODE2_A2_ROW);
xlabel('$a_5$','Interpreter','latex','FontSize',16);
ylabel('$a_2$','Interpreter','latex','FontSize',16);
title('ODE','FontSize',15,'FontWeight','bold'); xlim(lim.a52); ylim(lim.a2); grid on; pbaspect([1 1 1]);
addPanelLabel(figTP2, ax_t1, 'A');

ax_t2 = nexttile(tl6); hold on; box on; set(ax_t2,'Color','w','FontSize',12);
plotTwoParamLPA(lpa_hh, lpa_lp_segs, lim, LPA2_A5_ROW, LPA2_A2_ROW);
xlabel('$a_5$','Interpreter','latex','FontSize',16);
ylabel('$a_2$','Interpreter','latex','FontSize',16);
title('LPA','FontSize',15,'FontWeight','bold'); xlim(lim.a52); ylim(lim.a2); grid on; pbaspect([1 1 1]);
addPanelLabel(figTP2, ax_t2, 'B');

ax_t3 = nexttile(tl6); hold on; box on; set(ax_t3,'Color','w','FontSize',12);
plotTwoParamOverlay(ode_hh,ode_lp_segs,lpa_hh,lpa_lp_segs,lim,...
    ODE2_A5_ROW,ODE2_A2_ROW,LPA2_A5_ROW,LPA2_A2_ROW);
xlabel('$a_5$','Interpreter','latex','FontSize',16);
ylabel('$a_2$','Interpreter','latex','FontSize',16);
title('ODE vs LPA','FontSize',15,'FontWeight','bold'); xlim(lim.a52); ylim(lim.a2); grid on; pbaspect([1 1 1]);
addPanelLabel(figTP2, ax_t3, 'C');

exportgraphics(figTP2, fullfile(fig_folder,['TwoParam_ODE_LPA_Panels_abc',file_suffix,'.jpg']),'Resolution',300);

% =========================================================================
% 13. SIX-PANEL COMBINED FIGURE
% =========================================================================
figCombined = figure('Color','w','Units','inches',...
    'Position',[0 0 24 16],'PaperUnits','inches',...
    'PaperSize',[24 16],'PaperPosition',[0 0 24 16]);
tlC = tiledlayout(figCombined,2,3,'TileSpacing','compact','Padding','loose');

axC1 = nexttile(tlC); hold on; box on; set(axC1,'Color','w','FontSize',12);
plotBranch(lpa_ep,1,LPA1_A5_ROW,'b','r');
plotBranch(lpa_bp,1,LPA1_A5_ROW,'c','k');
plotLC_Envelope(lpa_lc_cells,1,4,'g');
labelBifurcations(lpa_ep,1,LPA1_A5_ROW);
xlabel('$a_5$','Interpreter','latex','FontSize',16);
ylabel('$R$','Interpreter','latex','FontSize',16);
xlim(lim.a5); ylim(lim.y_lpa); grid on; pbaspect([1 1 1]);
addPanelLabel(figCombined, axC1, 'A');

axC2 = nexttile(tlC); hold on; box on; set(axC2,'Color','w','FontSize',12);
plotBranch(lpa_ep,2,LPA1_A5_ROW,'b','r');
plotBranch(lpa_bp,2,LPA1_A5_ROW,'c','k');
plotLC_Envelope(lpa_lc_cells,2,4,'g');
labelBifurcations(lpa_ep,2,LPA1_A5_ROW);
xlabel('$a_5$','Interpreter','latex','FontSize',16);
ylabel('$\rho$','Interpreter','latex','FontSize',16);
xlim(lim.a5); ylim(lim.y_lpa); grid on; pbaspect([1 1 1]);
addPanelLabel(figCombined, axC2, 'B');

axC3 = nexttile(tlC); hold on; box on; set(axC3,'Color','w','FontSize',12);
plotTwoParamOverlay(ode_hh,ode_lp_segs,lpa_hh,lpa_lp_segs,lim,...
    ODE2_A5_ROW,ODE2_A2_ROW,LPA2_A5_ROW,LPA2_A2_ROW);
xlabel('$a_5$','Interpreter','latex','FontSize',16);
ylabel('$a_2$','Interpreter','latex','FontSize',16);
xlim(lim.a52); ylim(lim.a2); grid on; pbaspect([1 1 1]);
addPanelLabel(figCombined, axC3, 'C');

axC4 = nexttile(tlC); hold on; box on; set(axC4,'Color','w','FontSize',12);
plotBranch(ode_ep,1,ODE1_A5_ROW,'b','r');
plotLC_Envelope(ode_lc_cells,1,2,'g');
labelBifurcations(ode_ep,1,ODE1_A5_ROW);
xlabel('$a_5$','Interpreter','latex','FontSize',16);
ylabel('$R$','Interpreter','latex','FontSize',16);
xlim(lim.a5); ylim(lim.y_ode); grid on; pbaspect([1 1 1]);
addPanelLabel(figCombined, axC4, 'D');

axC5 = nexttile(tlC); hold on; box on; set(axC5,'Color','w','FontSize',12);
plotBranch(ode_ep,2,ODE1_A5_ROW,'b','r');
plotLC_Envelope(ode_lc_cells,2,2,'g');
labelBifurcations(ode_ep,2,ODE1_A5_ROW);
xlabel('$a_5$','Interpreter','latex','FontSize',16);
ylabel('$\rho$','Interpreter','latex','FontSize',16);
xlim(lim.a5); ylim(lim.y_ode); grid on; pbaspect([1 1 1]);
addPanelLabel(figCombined, axC5, 'E');

axC6 = nexttile(tlC); hold on; box on; set(axC6,'Color','w','FontSize',12);
plotTwoParamODE(ode_hh, ode_lp_segs, lim, ODE2_A5_ROW, ODE2_A2_ROW);
xlabel('$a_5$','Interpreter','latex','FontSize',16);
ylabel('$a_2$','Interpreter','latex','FontSize',16);
xlim(lim.a52); ylim(lim.a2); grid on; pbaspect([1 1 1]);
addPanelLabel(figCombined, axC6, 'F');

exportgraphics(figCombined, ...
    fullfile(fig_folder,['Combined_LPA_GB_SixPanel',file_suffix,'.jpg']), ...
    'Resolution',300);

fprintf('BP_BP figures skipped — BP_BP data not loaded.\n');

% =========================================================================
% HELPER FUNCTIONS
% =========================================================================

function data = combineHH(path, base, a5r, a2r)
% Loads all H_H(k).mat files in consecutive pairs and stitches into one
% ordered curve. MatCont saves each continuation restart as a new pair:
% (1)+(2) first restart, (3)+(4) second, (5)+(6) third, etc.
% Within each pair, file (k+1) is the backward direction and is flipped
% before prepending to file (k), exactly as combineFiles does.
% All pairs are then concatenated in order.
    data.x = []; data.f = []; data.s = [];
    k = 1;
    while true
        f1 = fullfile(path, sprintf('%s(%d).mat', base, k));
        f2 = fullfile(path, sprintf('%s(%d).mat', base, k+1));
        if ~exist(f1, 'file'), break; end

        d1 = load(f1);
        seg.x = d1.x; seg.f = []; seg.s = [];
        if isfield(d1,'f'), seg.f = d1.f; end
        if isfield(d1,'s'), seg.s = d1.s; end

        if exist(f2, 'file')
            d2 = load(f2);
            % Flip file (k+1) and prepend to file (k)
            x2 = fliplr(d2.x); x2(:,end) = [];
            seg.x = [x2, d1.x];
            if isfield(d1,'f') && isfield(d2,'f')
                f2x = fliplr(d2.f); f2x(:,end) = [];
                seg.f = [f2x, d1.f];
            end
            % Merge special points from both halves
            off = size(x2, 2);
            s1 = []; s2 = [];
            if isfield(d1,'s') && ~isempty(d1.s)
                s1 = d1.s;
                for j = 1:length(s1), s1(j).index = s1(j).index + off; end
            end
            if isfield(d2,'s') && ~isempty(d2.s)
                N = size(d2.x, 2);
                for j = 1:length(d2.s), d2.s(j).index = N - d2.s(j).index + 1; end
                s2 = d2.s([d2.s.index] < N);
            end
            seg.s = [s2(:); s1(:)];
        end

        % Concatenate this pair onto the growing curve
        if isempty(data.x)
            data.x = seg.x;
            data.f = seg.f;
            data.s = seg.s(:);
        else
            off2 = size(data.x, 2);
            data.x = [data.x, seg.x];
            if ~isempty(seg.f), data.f = [data.f, seg.f]; end
            if ~isempty(seg.s)
                ts = seg.s;
                for j = 1:length(ts), ts(j).index = ts(j).index + off2; end
                data.s = [data.s(:); ts(:)];
            end
        end

        k = k + 2;
    end
    fprintf('  combineHH: %d file(s) -> %d pts for %s\n', k-1, size(data.x,2), base);
end

function c = blendWithWhite(color, alpha)
    c = alpha * color(:)' + (1 - alpha) * [1, 1, 1];
    c = min(max(c, 0), 1);
end

function fillFakeAlpha(xA, yA, cA, xB, yB, cB, alpha)
    cAb   = blendWithWhite(cA, alpha);
    cBb   = blendWithWhite(cB, alpha);
    cOvlp = blendWithWhite(min(cA + cB, 1), alpha);
    fill(xA,yA,cAb,'EdgeColor','none','FaceAlpha',1,'HandleVisibility','off');
    fill(xB,yB,cBb,'EdgeColor','none','FaceAlpha',1,'HandleVisibility','off');
    try
        psA = polyshape(xA(:),yA(:),'Simplify',true,'KeepCollinearPoints',false);
        psB = polyshape(xB(:),yB(:),'Simplify',true,'KeepCollinearPoints',false);
        psI = intersect(psA,psB);
        if psI.NumRegions > 0
            [xi,yi] = boundary(psI);
            fill(xi,yi,cOvlp,'EdgeColor','none','FaceAlpha',1,'HandleVisibility','off');
        end
    catch
        fill(xB,yB,cB,'EdgeColor','none','FaceAlpha',alpha,'HandleVisibility','off');
    end
end

function data = combineFiles(path, base)
    data.x=[]; data.f=[]; data.s=[];
    f1=fullfile(path,sprintf('%s(1).mat',base));
    f2=fullfile(path,sprintf('%s(2).mat',base));
    h1=exist(f1,'file'); h2=exist(f2,'file');
    if ~h1&&~h2, return; end
    if h1&&~h2
        tmp=load(f1); data.x=tmp.x;
        if isfield(tmp,'f'),data.f=tmp.f;end
        if isfield(tmp,'s'),data.s=tmp.s;end
    else
        d1=load(f1); d2=load(f2);
        x2=fliplr(d2.x); x2(:,end)=[];
        f2x=[]; if isfield(d2,'f'),f2x=fliplr(d2.f);f2x(:,end)=[];end
        s2=[];
        if isfield(d2,'s')&&~isempty(d2.s)
            N=size(d2.x,2);
            for j=1:length(d2.s), d2.s(j).index=N-d2.s(j).index+1;end
            s2=d2.s([d2.s.index]<N);
        end
        data.x=[x2,d1.x];
        if ~isempty(f2x)&&isfield(d1,'f'),data.f=[f2x,d1.f];end
        off=size(x2,2); s1=[];if isfield(d1,'s'),s1=d1.s;end
        for j=1:length(s1),s1(j).index=s1(j).index+off;end
        data.s=[s2(:);s1(:)];
    end
end

function data = combineAllFiles(path, base)
    data.x=[]; data.f=[]; data.s=[]; i=1; off=0;
    while true
        fn=fullfile(path,sprintf('%s(%d).mat',base,i));
        if ~exist(fn,'file'),break;end
        tmp=load(fn); data.x=[data.x,tmp.x];
        if isfield(tmp,'f'),data.f=[data.f,tmp.f];end
        if isfield(tmp,'s')&&~isempty(tmp.s)
            ts=tmp.s; for j=1:length(ts),ts(j).index=ts(j).index+off;end
            data.s=[data.s(:);ts(:)];
        end
        off=off+size(tmp.x,2); i=i+1;
    end
    if i>1, fprintf('  combineAllFiles: %d seg(s) for %s, %d pts\n',i-1,base,size(data.x,2));end
end

function cellData = loadIndependently(path, base)
    cellData={}; i=1;
    while true
        fn=fullfile(path,sprintf('%s(%d).mat',base,i));
        if ~exist(fn,'file'),break;end
        cellData{i}=load(fn); i=i+1;
    end
end

function plotBranch(data, varIdx, paramIdx, stableClr, unstableClr)
    if isempty(data.x),return;end
    xv=data.x(paramIdx,:); yv=data.x(varIdx,:);
    isu=max(real(data.f),[],1)>1e-8;
    for i=1:length(xv)-1
        if abs(xv(i+1)-xv(i))>2,continue;end
        if isu(i)||isu(i+1), clr=[unstableClr '--']; else, clr=[stableClr '-'];end
        plot(xv(i:i+1),yv(i:i+1),clr,'LineWidth',2);
    end
end

function plotLC_Envelope(cellData, varIdx, nDim, clr)
    if isempty(cellData),return;end
    for k=1:length(cellData)
        d=cellData{k}; xm=d.x; fm=d.f; a5=xm(end,:);
        mi=varIdx:nDim:size(xm,1)-2;
        lc_max=max(xm(mi,:),[],1); lc_min=min(xm(mi,:),[],1);
        isu=max(abs(fm(end-3:end,:)),[],1)>1.001;
        for i=1:length(a5)-1
            if isu(i)||isu(i+1),sty='m--';else,sty=[clr '-'];end
            plot(a5(i:i+1),lc_max(i:i+1),sty,'LineWidth',2);
            plot(a5(i:i+1),lc_min(i:i+1),sty,'LineWidth',2);
        end
        if isfield(d,'s')
            for j=1:length(d.s)
                lab = strtrim(d.s(j).label);
                if strcmp(lab,'NS')
                    idx=d.s(j).index;
                    line([a5(idx) a5(idx)],[lc_min(idx) lc_max(idx)],'Color','y','LineWidth',2.5);
                    text(a5(idx),lc_max(idx),' NS','FontSize',10,'FontWeight','bold','Color','y','VerticalAlignment','bottom');
                end
            end
        end
    end
end

function labelBifurcations(data, varIdx, paramIdx)
    if isempty(data.s),return;end
    bif=struct('x',{},'y',{},'lab',{},'mstyle',{},'lbl',{});
    for i=1:length(data.s)
        lab=strtrim(data.s(i).label);
        if ~ismember(lab,{'H','LP','BP','N'}),continue;end
        idx=data.s(i).index;
        if idx<1||idx>size(data.x,2),continue;end
        e.x=data.x(paramIdx,idx); e.y=data.x(varIdx,idx); e.lab=lab;
        if e.x < 0, continue; end
        switch lab
            case 'H',  e.mstyle='r*'; e.lbl='HB';
            case 'LP', e.mstyle='k*'; e.lbl='LP';
            case 'BP', e.mstyle='ro'; e.lbl='BP';
            case 'N',  e.mstyle='ms'; e.lbl='N';
        end
        bif(end+1)=e; %#ok<AGROW>
    end
    if isempty(bif),return;end
    xv=[bif.x]; yv=[bif.y];
    xr=range(xv); if xr<eps,xr=1;end
    yr=range(yv); if yr<eps,yr=0.1;end
    THR=0.15; va_d=zeros(1,numel(bif));
    for i=1:numel(bif)
        for j=i+1:numel(bif)
            if abs(bif(i).x-bif(j).x)/xr<THR&&abs(bif(i).y-bif(j).y)/yr<THR
                if bif(i).y>=bif(j).y,va_d(i)=1;va_d(j)=-1;
                else,va_d(i)=-1;va_d(j)=1;end
            end
        end
    end
    for i=1:numel(bif)
        plot(bif(i).x,bif(i).y,bif(i).mstyle,'MarkerSize',10,'LineWidth',2);
        if va_d(i)>=0,va='bottom';else,va='top';end
        text(bif(i).x,bif(i).y,[' ',bif(i).lbl],'FontSize',10,'FontWeight','bold',...
            'HorizontalAlignment','left','VerticalAlignment',va,'Color','k');
    end
end

function [x_out,y_out] = collectLPData(lpSegs, a5r, a2r)
    x_out=[]; y_out=[];
    if isempty(lpSegs),return;end
    s1=lpSegs{1};
    if ~isempty(s1.x), x_out=fliplr(s1.x(a5r,:)); y_out=fliplr(s1.x(a2r,:));end
    for k=2:length(lpSegs)
        sg=lpSegs{k}; if isempty(sg.x),continue;end
        x_out=[x_out,sg.x(a5r,2:end)]; %#ok<AGROW>
        y_out=[y_out,sg.x(a2r,2:end)]; %#ok<AGROW>
    end
    v=isfinite(x_out)&isfinite(y_out); x_out=x_out(v); y_out=y_out(v);
    [x_out,idx]=sort(x_out); y_out=y_out(idx);
end

function plotTwoParamODE(hh, lpSegs, lim, a5r, a2r)
    G=[0,1,0]; R=[1,0,0]; A=0.5;
    hasHH=~isempty(hh.x); hasLP=~isempty(lpSegs);
    xLP=[];yLP=[];xHH=[];yHH=[];
    if hasLP
        [xl,yl]=collectLPData(lpSegs,a5r,a2r);
        if numel(xl)>=2,xLP=[xl,fliplr(xl),xl(1)];yLP=[yl,fliplr(yl),yl(1)];end
    end
    if hasHH
        % No sort — preserve MatCont curve order
        x=hh.x(a5r,:); y=hh.x(a2r,:);
        v=isfinite(x)&isfinite(y); x=x(v); y=y(v);
        xHH=[x,fliplr(x),x(1)]; yHH=[y,fliplr(y),y(1)];
    end
    if ~isempty(xLP)&&~isempty(xHH),fillFakeAlpha(xLP,yLP,G,xHH,yHH,R,A);
    elseif ~isempty(xLP),fill(xLP,yLP,blendWithWhite(G,A),'EdgeColor','none','FaceAlpha',1,'HandleVisibility','off');
    elseif ~isempty(xHH),fill(xHH,yHH,blendWithWhite(R,A),'EdgeColor','none','FaceAlpha',1,'HandleVisibility','off');end
    if ~isempty(xLP)
        [xl,yl]=collectLPData(lpSegs,a5r,a2r);
        plot(xl,yl,'g-','LineWidth',2,'HandleVisibility','off');
    end
    if hasHH
        x=hh.x(a5r,:); y=hh.x(a2r,:);
        v=isfinite(x)&isfinite(y); x=x(v); y=y(v);
        plot(x,y,'r-','LineWidth',2,'HandleVisibility','off');
    end
end

function plotTwoParamLPA(hh, lpSegs, lim, a5r, a2r)
    G=[0,1,0]; R=[1,0,0]; A=0.5;
    hasHH=~isempty(hh.x); hasLP=~isempty(lpSegs);
    xLP=[];yLP=[];xHH=[];yHH=[];
    if hasLP
        [xl,yl]=collectLPData(lpSegs,a5r,a2r);
        if numel(xl)>=2,xLP=[xl,fliplr(xl),xl(1)];yLP=[yl,fliplr(yl),yl(1)];end
    end
    if hasHH
        % No sort — preserve MatCont curve order
        x=hh.x(a5r,:); y=hh.x(a2r,:);
        v=isfinite(x)&isfinite(y); x=x(v); y=y(v);
        xHH=[x,fliplr(x),x(1)]; yHH=[y,fliplr(y),y(1)];
    end
    if ~isempty(xLP)&&~isempty(xHH),fillFakeAlpha(xLP,yLP,G,xHH,yHH,R,A);
    elseif ~isempty(xLP),fill(xLP,yLP,blendWithWhite(G,A),'EdgeColor','none','FaceAlpha',1,'HandleVisibility','off');
    elseif ~isempty(xHH),fill(xHH,yHH,blendWithWhite(R,A),'EdgeColor','none','FaceAlpha',1,'HandleVisibility','off');end
    if ~isempty(xLP)
        [xl,yl]=collectLPData(lpSegs,a5r,a2r);
        plot(xl,yl,'g-','LineWidth',2,'HandleVisibility','off');
    end
    if hasHH
        x=hh.x(a5r,:); y=hh.x(a2r,:);
        v=isfinite(x)&isfinite(y); x=x(v); y=y(v);
        plot(x,y,'r-','LineWidth',2,'HandleVisibility','off');
    end
end

function plotTwoParamOverlay(ode_hh,ode_lpSegs,lpa_hh,lpa_lpSegs,lim,...
                              ode_a5r,ode_a2r,lpa_a5r,lpa_a2r)
    G=[0,1,0]; R=[1,0,0]; A=0.5;

    function [xP,yP,xR,yR] = buildHH(hh,a5r,a2r)
        xP=[]; yP=[]; xR=[]; yR=[];
        if isempty(hh.x), return; end
        % No sort — preserve MatCont curve order
        x=hh.x(a5r,:); y=hh.x(a2r,:);
        v=isfinite(x)&isfinite(y); x=x(v); y=y(v);
        xR=x; yR=y;
        xP=[x,fliplr(x),x(1)]; yP=[y,fliplr(y),y(1)];
    end

    xOLP=[];yOLP=[];xLLP=[];yLLP=[];
    if ~isempty(ode_lpSegs)
        [xl,yl]=collectLPData(ode_lpSegs,ode_a5r,ode_a2r);
        if numel(xl)>=2,xOLP=[xl,fliplr(xl),xl(1)];yOLP=[yl,fliplr(yl),yl(1)];end
    end
    if ~isempty(lpa_lpSegs)
        [xl,yl]=collectLPData(lpa_lpSegs,lpa_a5r,lpa_a2r);
        if numel(xl)>=2,xLLP=[xl,fliplr(xl),xl(1)];yLLP=[yl,fliplr(yl),yl(1)];end
    end
    [xOHH,yOHH,xOHHr,yOHHr]=buildHH(ode_hh,ode_a5r,ode_a2r);
    [xLHH,yLHH,xLHHr,yLHHr]=buildHH(lpa_hh,lpa_a5r,lpa_a2r);

    if ~isempty(xOLP)&&~isempty(xOHH),fillFakeAlpha(xOLP,yOLP,G,xOHH,yOHH,R,A);
    elseif ~isempty(xOLP),fill(xOLP,yOLP,blendWithWhite(G,A),'EdgeColor','none','FaceAlpha',1,'HandleVisibility','off');
    elseif ~isempty(xOHH),fill(xOHH,yOHH,blendWithWhite(R,A),'EdgeColor','none','FaceAlpha',1,'HandleVisibility','off');end

    if ~isempty(xLLP)&&~isempty(xLHH),fillFakeAlpha(xLLP,yLLP,G,xLHH,yLHH,R,A);
    elseif ~isempty(xLLP),fill(xLLP,yLLP,blendWithWhite(G,A),'EdgeColor','none','FaceAlpha',1,'HandleVisibility','off');
    elseif ~isempty(xLHH),fill(xLHH,yLHH,blendWithWhite(R,A),'EdgeColor','none','FaceAlpha',1,'HandleVisibility','off');end

    if ~isempty(ode_lpSegs)
        [xl,yl]=collectLPData(ode_lpSegs,ode_a5r,ode_a2r);
        if numel(xl)>=2,plot(xl,yl,'g-','LineWidth',2,'HandleVisibility','off');end
    end
    if ~isempty(lpa_lpSegs)
        [xl,yl]=collectLPData(lpa_lpSegs,lpa_a5r,lpa_a2r);
        if numel(xl)>=2,plot(xl,yl,'g--','LineWidth',2,'HandleVisibility','off');end
    end
    if ~isempty(xOHHr),plot(xOHHr,yOHHr,'r-','LineWidth',2,'HandleVisibility','off');end
    if ~isempty(xLHHr),plot(xLHHr,yLHHr,'r--','LineWidth',2,'HandleVisibility','off');end
end

function addPanelLabel(fig, ax, labelStr)
    drawnow('nocallbacks');
    text(ax, 0, 1.02, ['\bf' labelStr], ...
        'Units',               'normalized', ...
        'FontName',            'Times New Roman', ...
        'FontSize',            14, ...
        'FontWeight',          'bold', ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment',   'bottom', ...
        'Interpreter',         'tex', ...
        'Clipping',            'off');
end