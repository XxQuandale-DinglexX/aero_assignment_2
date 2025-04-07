clc;
clear;
%close all;
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
AR_values = [4, 6, 8, 10, 1e5];
alpha_range = linspace(-4, 10, 15)*pi/180;
alpha_lo = -4.203*pi/180;
beta = linspace(1e-5, pi, n);
x = (1 - cos(beta)) / 2;

data = load('NACA4415_Detailed.txt');
alpha_data = data(:,1)*pi/180;
Cl_data = data(:,2);
b_span = 1;
M0 = gradient(Cl_data, alpha_data);

%% Initialize Variables
A_n_alphas = zeros(length(alpha_range), n, length(AR_values));
Cl_values = zeros(length(alpha_range), length(AR_values));
CDi_values = zeros(length(alpha_range), length(AR_values));

%% Tapered Wing
TR_values = [0.2, 0.4, 0.6, 0.8, 1.0];
AR_tapered = 6;
alpha_tapered = 6;
%% Compute A_n, Cl, and CDi for Each Aspect Ratio
for i_AR = 1:length(AR_values)
    b_c = AR_values(i_AR);
    
    for i_alpha = 1:length(alpha_range)
        m0 = interp1(alpha_data, M0, alpha_lo, 'linear', 'extrap');
        alpha_i = alpha_range(i_alpha);
        
        coeffs = zeros(n, n);
        for i_theta = 1:n
            theta = beta(i_theta);
            C = -4 * b_c / m0;

            for i_an = 1:n
                coeffs(i_theta, i_an) = C * sin(i_an * theta) - (i_an * sin(i_an * theta) / sin(theta));
            end
        end
        
        % Solve for A_n coefficients
        b = (alpha_lo - alpha_i) * ones(n, 1);
        A = coeffs \ b;
        A_n_alphas(i_alpha, :, i_AR) = A';
    end
    
    % Compute Cl and CDi
    Cl_values(:, i_AR) = pi * b_c * A_n_alphas(:, 1, i_AR);
    
    delta = zeros(length(alpha_range), 1);
    for i = 2:n
        delta = delta + i * (A_n_alphas(:, i, i_AR) ./ A_n_alphas(:, 1, i_AR)).^2;
    end
    
    CDi_values(:, i_AR) = Cl_values(:, i_AR).^2 ./ (pi * b_c) .* (1 + delta);
end

%% Plot 1: C_L vs AoA
fontsize_general = 18;

figure(12), clf;
hold on;
xlim([-5, 10]);
ylim([-1, 2]);
plot(xlim, [0,0], 'LineWidth', 1, 'HandleVisibility', 'off', 'Color', 'k');
plot([0, 0], ylim, 'k', 'LineWidth', 1, 'HandleVisibility', 'off');

AR_values_str = {num2str(4), num2str(6), num2str(8), num2str(10), num2str(inf)};
colors = ["#00a5cf", "#d1495b", "#3c1642", "#4f772d", "#edae49"];
for i_AR = 1:length(AR_values)
    plot(alpha_range*180/pi, Cl_values(:, i_AR), 'LineWidth', 2, ...
        'DisplayName', ['AR = ' AR_values_str{i_AR}], 'Color', colors(i_AR));
end
grid on;
xticks(-4:2:10)

%title('Lift Coefficient vs Angle of Attack');
legend('Location', 'southeast', 'Interpreter', 'latex');

set(gca, 'Color', 'white', 'FontSize', fontsize_general ,'Xcolor', 'k', 'YColor', 'k', ...
    'GridColor', [0.1, 0.1, 0.1], 'XMinorGrid', 'off', 'YMinorGrid', 'off');

xlabel('AoA [deg]', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$C_L$', 'Interpreter', 'latex', 'FontSize', 22);
%% Plot 2: CDi vs AoA
figure(13), clf;
hold on;

xlim([-5, 10]);
ylim([-0.01, 0.1]);
xticks(-4:2:10)

plot(xlim, [0,0], 'LineWidth', 1, 'HandleVisibility', 'off', 'Color', 'k');
plot([0, 0], ylim, 'k', 'LineWidth', 1, 'HandleVisibility', 'off');
for i_AR = 1:length(AR_values)
    plot(alpha_range*180/pi, CDi_values(:, i_AR), 'LineWidth', 2, ...
        'DisplayName', ['AR = ' AR_values_str{i_AR}], 'Color', colors(i_AR));
end
grid on;

%title('Induced Drag Coefficient vs Angle of Attack');
legend('Location', 'northwest', 'Interpreter','latex');
set(gca, 'Color', 'white', 'FontSize', fontsize_general ,'Xcolor', 'k', 'YColor', 'k', ...
    'GridColor', [0.1, 0.1, 0.1], 'XMinorGrid', 'off', 'YMinorGrid', 'off');
xlabel('AoA [deg]', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$C_{D,i}$', 'Interpreter', 'latex', 'FontSize', 22);

%% Calculation for Induced AoA
target_alphas = [0, 5*pi/180, 10*pi/180];
alpha_i_values = zeros(length(AR_values), length(beta), length(target_alphas));

for i_target_alphas = 1:length(target_alphas)
    target_alpha = target_alphas(i_target_alphas);
    target_index = find(abs(alpha_range - target_alpha) < 1e-5, 1);

    for i_AR = 1:length(AR_values)
        A_n = A_n_alphas(target_index, :, i_AR); % Get A_n for specific AoA and AR
        
        for i_theta = 1:length(beta)
            theta = beta(i_theta);
            % Compute the induced angle of attack
            alpha_i_values(i_AR, i_theta, i_target_alphas) = sum((1:n) .* A_n .* sin((1:n) * theta) ./ sin(theta));
        end
    end
end

alpha_i_values = rad2deg(alpha_i_values);

%% Convert Theta to Spanwise Position y/b
y = (b_span / 2) * cos(beta);

%% Plot for each Target AoA (0°, 5°, 10°) with Different ARs
colors = ["#00a5cf", "#d1495b", "#3c1642", "#4f772d", "#edae49"];


for i_target_alphas = 1:length(target_alphas)
    target_alpha = target_alphas(i_target_alphas);

    figure(i_target_alphas), clf;
    hold on;
    
    
    xlim([-0.55, 0.55]);
    ylim([-0.02,15]);
    xticks(-0.5:0.1:0.5)
    yticks(0:2:14);
    
    %plot(xlim, [0,0], 'LineWidth', 1, 'HandleVisibility', 'off', 'Color', 'k');
    plot([0, 0], ylim, 'k', 'LineWidth', 1, 'HandleVisibility', 'off');

    for i_AR = 1:length(AR_values)
        plot(y, squeeze(alpha_i_values(i_AR, :, i_target_alphas)), 'LineWidth', 2, 'DisplayName', ['AR = ' AR_values_str{i_AR}], 'Color', colors(i_AR));
    end
    
    grid on;

    % title(['Induced Angle of Attack Along the Wingspan at \alpha = ' num2str(target_alpha*180/pi) '^\circ']);
    title(['AoA = ' num2str(target_alpha*180/pi) '^\circ']);

    legend('Location', 'best', 'Interpreter', 'latex');
    xlim([-0.55, 0.55]); % Restrict x-axis range

    set(gca, 'Color', 'white', 'FontSize', fontsize_general ,'Xcolor', 'k', 'YColor', 'k', ...
    'GridColor', [0.1, 0.1, 0.1], 'XMinorGrid', 'off', 'YMinorGrid', 'off');

    xlabel('$x/b$', 'Interpreter', 'latex', 'FontSize', 22);
    ylabel('$\alpha_i$ [deg]', 'Interpreter', 'latex', 'FontSize', 22);
end