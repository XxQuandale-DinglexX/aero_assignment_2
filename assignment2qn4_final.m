clc;
clear;
close all;

%% Set default figure properties to light mode
close all;
set(0, 'DefaultFigureColor', 'w');     % Default figure background color (white)
set(0, 'DefaultAxesColor', 'w');   % Transparent axes background
set(0, 'DefaultAxesXColor', 'k');     % Default black X-axis
set(0, 'DefaultAxesYColor', 'k');     % Default black Y-axis
set(0, 'DefaultAxesGridColor', [0.15, 0.15, 0.15]); % Default grid line color
set(0, 'DefaultTextColor', 'k');      % Default text color (black)
set(0, 'DefaultFigurePosition',  [80, 50, 750, 450]);

%% Airfoil Parameters
n = 500;
b_c = 6; % aspect ratio
b_span = 1;

alpha_g0 = 4 * pi / 180;
Alpha_tipvalues = [0, 2, 4, 6, 8] * pi / 180; 

alpha_lo = -4.203 * pi / 180;
beta = linspace(1e-5, pi, n);
x = (1 - cos(beta)) / 2;

% Load airfoil data
data = load('NACA4415_Detailed.txt');
alpha_data = data(:,1) * pi / 180;
Cl_data = data(:,2);
M0 = gradient(Cl_data, alpha_data);
Q_inf = 1;

% Initialize variables
A_n_alphas = zeros(n, n, length(Alpha_tipvalues));
Cl_values = zeros(n, length(Alpha_tipvalues));
CDi_values = zeros(n, length(Alpha_tipvalues));
alpha_i_values = zeros(length(Alpha_tipvalues), n);
Gamma_values = zeros(length(Alpha_tipvalues), n);

% Loop over geometric tip angles
for i_tip = 1:length(Alpha_tipvalues)
    alpha_g_tip = Alpha_tipvalues(i_tip);
    m0 = interp1(alpha_data, M0, alpha_lo, 'linear', 'extrap');
    C = -4 * b_c / m0;

    % Compute geometric AoA distribution across span
    x_bar = -cos(beta);
    Geom_alpha_range = alpha_g0 + (alpha_g_tip - alpha_g0) * abs(x_bar);

    coeffs = zeros(n, n);

    for i_theta = 1:n
        theta = beta(i_theta);
        for i_an = 1:n
            coeffs(i_theta, i_an) = C * sin(i_an * theta) - (i_an * sin(i_an * theta) / sin(theta));
        end
    end

    % Solve for A_n
    for i_alpha = 1:n
        b = alpha_lo - Geom_alpha_range(:);
        A = coeffs \ b;
        A_n_alphas(i_alpha,:,i_tip) = A';
    end

    % Compute Cl, CDi, alpha_i, and Gamma
    for i_alpha = 1:n
        A_n = A_n_alphas(i_alpha,:,i_tip);

        % Cl
        Cl_values(i_alpha, i_tip) = 4 * b_c * sum(A_n .* sin((1:n) .* beta(i_alpha)));

        % Induced AoA
        alpha_i_values(i_tip, i_alpha) = sum((1:n) .* A_n .* sin((1:n) .* beta(i_alpha)) ./ sin(beta(i_alpha)));

        % CDi
        CDi_values(i_alpha, i_tip) = Cl_values(i_alpha, i_tip) * alpha_i_values(i_tip, i_alpha);

        % Gamma
        Gamma_values(i_tip, i_alpha) = 2 * b_c * sum(A_n .* sin((1:n) * beta(i_alpha)));
    end
end

alpha_i_values = rad2deg(alpha_i_values);

% Convert theta to y/b
% y = (b_span / 2) * cos(beta);
alpha_geom_values_str = {'0', '2', '4', '6', '8'};
colors = ["#00a5cf", "#d1495b", "#3c1642", "#4f772d", "#edae49"];

%% Plot Cl vs \tilde{x}
fontsize_general = 18;
legendsize = 19;
figure(1), clf;
hold on;
for i = 1:length(Alpha_tipvalues)
    plot(x_bar, Cl_values(:, i), 'LineWidth', 2, ...
        'DisplayName', ['$\alpha_{g,tip}$ = ' alpha_geom_values_str{i} '$^{\circ}$' ], ...
        'Color', colors(i));
end
grid on;
set(gca, 'Color', 'white', 'FontSize', fontsize_general ,'Xcolor', 'k', 'YColor', 'k', ...
    'GridColor', [0.1, 0.1, 0.1], 'XMinorGrid', 'off', 'YMinorGrid', 'off');
xlim([-1, 1]);
xticks(-1:0.2:1);
plot(xlim, [0,0], 'LineWidth', 1, 'HandleVisibility', 'off', 'Color', 'k');
plot([0, 0], ylim, 'k', 'LineWidth', 1, 'HandleVisibility', 'off');



ylabel('$C_l$', 'Interpreter', 'latex', 'FontSize', 22);
xlabel('$\tilde{x}$', 'Interpreter', 'latex', 'FontSize', 22);

legend('Location', 'south', 'Interpreter', 'latex', 'FontSize', legendsize);


%%
% todo
% font size bigger
% legend with degrees

% task 2:
% remove ylabel
% adapt xlabel
%% Plot CDi vs \tilde{x}
figure(2), clf;
hold on;
for i = 1:length(Alpha_tipvalues)
    plot(x_bar, CDi_values(:, i), 'LineWidth', 2, ...
        'DisplayName', ['$\alpha_{g,tip}$ = ' alpha_geom_values_str{i} '$^{\circ}$'], ...
        'Color', colors(i));
end
grid on;
xlabel('$\tilde{x}$', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$C_{d,i}$', 'Interpreter', 'latex', 'FontSize', 22);
set(gca, 'Color', 'white', 'FontSize', fontsize_general ,'Xcolor', 'k', 'YColor', 'k', ...
    'GridColor', [0.1, 0.1, 0.1], 'XMinorGrid', 'off', 'YMinorGrid', 'off');
xlim([-1, 1]);
ylim([0, 0.075])
xticks(-1:0.2:1);
plot(xlim, [0,0], 'LineWidth', 1, 'HandleVisibility', 'off', 'Color', 'k');
plot([0, 0], ylim, 'k', 'LineWidth', 1, 'HandleVisibility', 'off');

legend('Location', 'north', 'Interpreter','latex', 'FontSize', legendsize);

%% Plot Induced AoA vs \tilde{x}
figure(3), clf;
hold on;
for i = 1:length(Alpha_tipvalues)
    plot(x_bar, alpha_i_values(i,:), 'LineWidth', 2, ...
        'DisplayName', ['$\alpha_{g,tip}$ = ' alpha_geom_values_str{i} '$^{\circ}$'], ...
        'Color', colors(i));
end
grid on;
xlabel('$\tilde{x}$', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$\alpha_i$ [deg]', 'Interpreter', 'latex', 'FontSize', 22);

set(gca, 'Color', 'white', 'FontSize', fontsize_general ,'Xcolor', 'k', 'YColor', 'k', ...
    'GridColor', [0.1, 0.1, 0.1], 'XMinorGrid', 'off', 'YMinorGrid', 'off');
xlim([-1, 1]);
xticks(-1:0.2:1);
plot(xlim, [0,0], 'LineWidth', 1, 'HandleVisibility', 'off', 'Color', 'k');
plot([0, 0], ylim, 'k', 'LineWidth', 1, 'HandleVisibility', 'off');



legend('Location', 'north', 'Interpreter','latex', 'FontSize', legendsize);
%% Plot Gamma vs \tilde{x}
figure(4), clf;
hold on;
for i = 1:length(Alpha_tipvalues)
    plot(x_bar, Gamma_values(i,:), 'LineWidth', 2, ...
        'DisplayName', ['$\alpha_{g,tip}$ = ' alpha_geom_values_str{i} '$^{\circ}$'], ...
        'Color', colors(i));
end
grid on;
xlabel('$\tilde{x}$', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$\tilde{\Gamma}$', 'Interpreter', 'latex', 'FontSize', 22);

set(gca, 'Color', 'white', 'FontSize', fontsize_general ,'Xcolor', 'k', 'YColor', 'k', ...
    'GridColor', [0.1, 0.1, 0.1], 'XMinorGrid', 'off', 'YMinorGrid', 'off');
xlim([-1, 1]);
xticks(-1:0.2:1);
plot(xlim, [0,0], 'LineWidth', 1, 'HandleVisibility', 'off', 'Color', 'k');
plot([0, 0], ylim, 'k', 'LineWidth', 1, 'HandleVisibility', 'off');


legend('Location', 'south', 'Interpreter','latex', 'FontSize', legendsize);