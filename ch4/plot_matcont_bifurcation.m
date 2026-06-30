%% Plot MatCont bifurcation output
clc; clear; close all

ep1 = load("EP_EP(1).mat");
ep2 = load("EP_EP(2).mat");

% stitch the two continuation branches together
R1_prep  = fliplr(ep1.x(1, 2:end));
rho1_prep  = fliplr(ep1.x(2, 2:end));
a51_prep = fliplr(ep1.x(3, 2:end));

R2  = ep2.x(1, :);
rho2  = ep2.x(2, :);
a52 = ep2.x(3, :);

R  = [R1_prep, R2];
rho  = [rho1_prep, rho2];
a5 = [a51_prep, a52];

f1_prep = fliplr(ep1.f(:, 2:end));
f2      = ep2.f;
f_combined = [f1_prep, f2];
re_f = real(f_combined);

index_shift = length(R1_prep);

% find limit point / Hopf indices flagged by MatCont
ep1_ncols = size(ep1.x, 2);

tp_msg = {}; tp_idx = [];
for i = 1:length(ep1.s)
    msg = lower(char(ep1.s(i).msg));
    orig_idx = double(ep1.s(i).index);
    if orig_idx < 2, continue; end
    if contains(msg, 'limit point') || contains(msg, 'hopf')
        comb_idx = ep1_ncols - orig_idx + 1;
        tp_msg{end+1} = msg;
        tp_idx(end+1) = comb_idx;
    end
end

h_indices = tp_idx(contains(tp_msg, 'hopf'));
idx_H1 = min(h_indices);
idx_H2 = max(h_indices);

seg1 = 1:idx_H1;
seg2 = idx_H1:idx_H2;
seg3 = idx_H2:length(a5);

%% Figure 1: one-parameter bifurcation diagrams
fig1 = figure('Name','Bifurcation Diagrams','Color','w','Position',[100 100 1200 500]);

subplot(1,2,1); hold on; box on;
plot(a5(seg1), R(seg1), 'b-', 'LineWidth', 2);
plot(a5(seg2), R(seg2), 'r--', 'LineWidth', 2);
plot(a5(seg3), R(seg3), 'b-', 'LineWidth', 2);
xlabel('a5'); ylabel('R');
set(gca,'FontSize',20);
xlim([0,4]); ylim([0,1]);

subplot(1,2,2); hold on; box on;
plot(a5(seg1), rho(seg1), 'b-', 'LineWidth', 2);
plot(a5(seg2), rho(seg2), 'r--', 'LineWidth', 2);
plot(a5(seg3), rho(seg3), 'b-', 'LineWidth', 2);
xlabel('a5'); ylabel('\rho');
set(gca,'FontSize',20);
xlim([0,4]); ylim([0,1]);

% mark Hopf (H) and limit points (LP)
for i = 1:length(tp_msg)
    idx = tp_idx(i);
    p_a5 = a5(idx);
    p_R = R(idx);
    p_rho = rho(idx);

    if contains(tp_msg{i}, 'hopf')
        label = 'H'; style = 'r*'; sizeM = 10;
    else
        label = 'LP'; style = 'k.'; sizeM = 25;
    end

    subplot(1,2,1);
    plot(p_a5, p_R, style, 'MarkerSize', sizeM);
    text(p_a5, p_R, ['  ', label], 'FontSize', 14);

    subplot(1,2,2);
    plot(p_a5, p_rho, style, 'MarkerSize', sizeM);
    text(p_a5, p_rho, ['  ', label], 'FontSize', 14);
end

exportgraphics(fig1,'Bifurcation_Diagrams.png','Resolution',600);

%% Two-parameter Hopf continuation (a5 vs a2)
h1_c1 = load("H_H(1)_c1_a2.mat");
h2_c1 = load("H_H(2)_c1_a2.mat");
h1_c2 = load("H_H(1)_c2_a2.mat");
h2_c2 = load("H_H(2)_c2_a2.mat");

sp_gh = []; sp_bt = [];
for i = 1:length(h1_c1.s)
    msg = lower(char(h1_c1.s(i).msg));
    idx = h1_c1.s(i).index;
    if contains(msg,'generalized hopf')
        sp_gh = [sp_gh; h1_c1.x(4,idx), h1_c1.x(3,idx)];
    elseif contains(msg,'bogdanov')
        sp_bt = [sp_bt; h1_c1.x(4,idx), h1_c1.x(3,idx)];
    end
end

a2_c1 = [fliplr(h1_c1.x(3,2:end)), h2_c1.x(3,:)];
a5_c1 = [fliplr(h1_c1.x(4,2:end)), h2_c1.x(4,:)];

a2_c2 = [fliplr(h1_c2.x(3,2:end)), h2_c2.x(3,:)];
a5_c2 = [fliplr(h1_c2.x(4,2:end)), h2_c2.x(4,:)];

fig2 = figure('Name','Two-Parameter Hopf','Color','w','Position',[150 150 800 600]);
hold on; box on;

valid_c1 = (a5_c1>=0)&(a5_c1<=40);
valid_c2 = (a5_c2>=0)&(a5_c2<=40);

a5_c1_f = a5_c1(valid_c1); a2_c1_f = a2_c1(valid_c1);
a5_c2_f = a5_c2(valid_c2); a2_c2_f = a2_c2(valid_c2);

[a5_c1_f,i1] = sort(a5_c1_f); a2_c1_f = a2_c1_f(i1);
[a5_c2_f,i2] = sort(a5_c2_f); a2_c2_f = a2_c2_f(i2);

% shade the oscillatory region between the two Hopf curves
X_fill = [a5_c1_f, fliplr(a5_c2_f)];
Y_fill = [a2_c1_f, fliplr(a2_c2_f)];

hOsc = fill(X_fill, Y_fill, 'r', 'EdgeColor','none', 'FaceAlpha',0.85, ...
            'DisplayName','Oscillatory');

hStable = fill(nan,nan,'w','EdgeColor','k','DisplayName','Stable');

plot(a5_c1, a2_c1, 'b-', 'LineWidth',2);
plot(a5_c2, a2_c2, 'b-', 'LineWidth',2);

xlim([0 17]);
ylim([2.5 15]);
xlabel('a_5');
ylabel('a_2');
set(gca,'FontSize',20);
legend([hOsc,hStable],'Location','northeast','FontSize',14);

exportgraphics(fig2,'Two_Param_Hopf_Bifurcation.png','Resolution',600);
