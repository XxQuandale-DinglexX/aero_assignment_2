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
fontsize_general = 18;

%% Airfoil Parameters
n = 500;
alpha_lo = -4.203*pi/180;
beta = linspace(1e-5, pi-1e-5, n);
x = (1 - cos(beta)) / 2;
alpha_i = deg2rad(6);

data = load('NACA4415_Detailed.txt');
alpha_data = data(:,1)*pi/180;
Cl_data = data(:,2);
M0 = gradient(Cl_data, alpha_data);
m0 = interp1(alpha_data, M0, alpha_lo, 'linear', 'extrap');

%% Tapered Wing
TR_values = [0.2, 0.4, 0.6, 0.8, 1.0];
AR_tapered = 6;
alpha_tapered = 6;
Cr = 1;

%% Initialize Variables
A_n_alphas = zeros(n, length(TR_values));

%% Compute A_n, Cl, and CDi for Each Aspect Ratio
for i_TR = 1:length(TR_values)
    TR = TR_values(i_TR);

    b_c = AR_tapered/2*Cr*(TR+1);

    coeffs = zeros(n, n);

    for i_theta = 1:n
        theta = beta(i_theta);
        x_bar = (- cos(theta));
        chord = Cr*((TR-1)*abs(x_bar)+1);
        C = -4 * b_c / m0 / chord;
        for i_an = 1:n
            coeffs(i_theta, i_an) = C * sin(i_an * theta) - (i_an * sin(i_an * theta) / sin(theta));
        end
    end
    
    % Solve for A_n coefficients
    b = (alpha_lo - alpha_i) * ones(n, 1);
    A = coeffs \ b;
    A_n_alphas(:, i_TR) = A';  % n, TR

end

%% Calculation for dimensionless Gamma_bar

Gamma_bar = zeros(length(beta), length(i_TR));
alpha_i_values = zeros(length(beta), length(i_TR));
Cl_values = zeros(length(beta), length(i_TR));
Cdi_values = zeros(length(beta), length(i_TR));

x_bar_list = -cos(beta);

for i_TR = 1:length(TR_values)
    TR = TR_values(i_TR);
    c_bar = Cr*(TR+1)/2;  
    A_n = A_n_alphas(:, i_TR)';
        
    for i_theta = 1:length(beta)
        theta = beta(i_theta);
        x_bar = (- cos(theta));
        chord = Cr*((TR-1)*abs(x_bar)+1);

        % Compute the induced angle of attack
        Gamma_bar(i_theta, i_TR) = AR_tapered*Cr*(TR+1)/c_bar*sum(A_n .* sin((1:n) * theta));
        alpha_i_values(i_theta, i_TR) = sum((1:n) .* A_n .* sin((1:n) * theta) ./ sin(theta));
        Cl_values(i_theta, i_TR) = 2*AR_tapered*Cr*(TR+1)/chord * sum(A_n .* sin((1:n) * theta));
        temp_for_Cdi = sum((1:n) .* A_n .* sin((1:n)*theta)/sin(theta));
        Cdi_values(i_theta, i_TR) = 2*AR_tapered*Cr*(TR+1)/chord * sum(A_n .* sin((1:n)*theta) * temp_for_Cdi);
    end
end

%% Plot Gamma
colors = ["#00a5cf", "#d1495b", "#3c1642", "#4f772d", "#edae49"];

figure(1), clf;
hold on;
for i_TR = 1:length(TR_values)
    plot(x_bar_list, Gamma_bar(:, i_TR), 'LineWidth', 2, ...
        'DisplayName', ['TR = ', num2str(TR_values(i_TR))], 'Color', colors(i_TR));
end
grid on;

legend('Location', 'south', 'Interpreter', 'latex');
set(gca, 'Color', 'white', 'FontSize', fontsize_general ,'Xcolor', 'k', 'YColor', 'k', ...
    'GridColor', [0.1, 0.1, 0.1], 'XMinorGrid', 'off', 'YMinorGrid', 'off');
xlim([-1, 1]);
xticks(-1:0.2:1);
plot(xlim, [0,0], 'LineWidth', 1, 'HandleVisibility', 'off', 'Color', 'k');
plot([0, 0], ylim, 'k', 'LineWidth', 1, 'HandleVisibility', 'off');


xlabel('$\tilde{x}$', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$\tilde{\Gamma}$', 'Interpreter', 'latex', 'FontSize', 22);
%% Plot ai
figure(2), clf;
hold on;
for i_TR = 1:length(TR_values)
    plot(x_bar_list, rad2deg(alpha_i_values(:, i_TR)), 'LineWidth', 2, ...
        'DisplayName', ['TR = ', num2str(TR_values(i_TR))], 'Color', colors(i_TR));
end
grid on;

legend('Location', 'best', 'Interpreter', 'latex');
set(gca, 'Color', 'white', 'FontSize', fontsize_general ,'Xcolor', 'k', 'YColor', 'k', ...
    'GridColor', [0.1, 0.1, 0.1], 'XMinorGrid', 'off', 'YMinorGrid', 'off');
xlim([-1, 1]);
xticks(-1:0.2:1);
plot(xlim, [0,0], 'LineWidth', 1, 'HandleVisibility', 'off', 'Color', 'k');
plot([0, 0], ylim, 'k', 'LineWidth', 1, 'HandleVisibility', 'off');


xlabel('$\tilde{x}$', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$\alpha_i$ [deg]', 'Interpreter', 'latex', 'FontSize', 22);
%% Plot cl
figure(3), clf;
hold on;
for i_TR = 1:length(TR_values)
    plot(x_bar_list, Cl_values(:, i_TR), 'LineWidth', 2, ...
        'DisplayName', ['TR = ', num2str(TR_values(i_TR))], 'Color', colors(i_TR));
end
grid on;

legend('Location', 'south', 'Interpreter', 'latex');
set(gca, 'Color', 'white', 'FontSize', fontsize_general ,'Xcolor', 'k', 'YColor', 'k', ...
    'GridColor', [0.1, 0.1, 0.1], 'XMinorGrid', 'off', 'YMinorGrid', 'off');
xlim([-1, 1]);
xticks(-1:0.2:1);
plot(xlim, [0,0], 'LineWidth', 1, 'HandleVisibility', 'off', 'Color', 'k');
plot([0, 0], ylim, 'k', 'LineWidth', 1, 'HandleVisibility', 'off');


ylabel('$C_l$', 'Interpreter', 'latex', 'FontSize', 22);
xlabel('$\tilde{x}$', 'Interpreter', 'latex', 'FontSize', 22);
%% Plot cdi
figure(4), clf;
hold on;
for i_TR = 1:length(TR_values)
    plot(x_bar_list, Cdi_values(:, i_TR), 'LineWidth', 2, ...
        'DisplayName', ['TR = ', num2str(TR_values(i_TR))], 'Color', colors(i_TR));
end
grid on;


legend('Location', 'south', 'Interpreter','latex');
set(gca, 'Color', 'white', 'FontSize', fontsize_general ,'Xcolor', 'k', 'YColor', 'k', ...
    'GridColor', [0.1, 0.1, 0.1], 'XMinorGrid', 'off', 'YMinorGrid', 'off');
xlim([-1, 1]);
xticks(-1:0.2:1);
plot(xlim, [0,0], 'LineWidth', 1, 'HandleVisibility', 'off', 'Color', 'k');
plot([0, 0], ylim, 'k', 'LineWidth', 1, 'HandleVisibility', 'off');


xlabel('$\tilde{x}$', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$C_{d,i}$', 'Interpreter', 'latex', 'FontSize', 22);
