clear all;
close all;

%% Airfoil Parameters

AR_values = [4, 6, 8, 10, inf];

data = load('NACA4415_2d_airfoil.txt');
alpha_data_complete = data(:,1);
Cl_data_complete = data(:,2);
Cd_data_complete = data(:,3);

target_index_list = find(-6<=alpha_data_complete & alpha_data_complete <=10);
alpha_lo = interp1(Cl_data_complete, alpha_data_complete, 0);

alpha_data = alpha_data_complete(target_index_list);
Cl_data = Cl_data_complete(target_index_list);
%% Initialize Variables

Cl_values = zeros(length(alpha_data), length(AR_values));
CDi_values = zeros(length(alpha_data), length(AR_values));
ai_values = zeros(length(alpha_data), length(AR_values));
CDf_values = zeros(length(alpha_data), length(AR_values));

%% Compute ai, CDi, CDf for Each Aspect Ratio

for i_AR = 1:length(AR_values)
    AR = AR_values(i_AR);
    Cl_values(:,i_AR) = 2*pi/(1+2/AR).*(deg2rad(alpha_data)-deg2rad(alpha_lo));
    CDi_values(:,i_AR) = 1/(pi*AR).*Cl_values(:,i_AR).^2;

    ai_values(:,i_AR) = rad2deg(2/(AR+2)*(deg2rad(alpha_data)-deg2rad(alpha_lo)));
    
    alpha_eff = 1/(1+2/AR).*alpha_data + 2/(AR+2)*alpha_lo;
    CDf_values(:,i_AR) = interp1(alpha_data_complete, Cd_data_complete, alpha_eff);
end


%% Set default figure properties to light mode
close all;
set(0, 'DefaultFigureColor', 'w');     % Default figure background color (white)
set(0, 'DefaultAxesColor', 'w');   % Transparent axes background
set(0, 'DefaultAxesXColor', 'k');     % Default black X-axis
set(0, 'DefaultAxesYColor', 'k');     % Default black Y-axis
set(0, 'DefaultAxesGridColor', [0.15, 0.15, 0.15]); % Default grid line color
set(0, 'DefaultTextColor', 'k');      % Default text color (black)
set(0, 'DefaultFigurePosition',  [80, 50, 750, 450]);

%% Plot 1: Cl vs AoA
% figure(1), clf;
% hold on;
% colors = ["#00a5cf", "#d1495b", "#3c1642", "#4f772d", "#edae49"];
% for i_AR = 1:length(AR_values)
%     plot(alpha_data, Cl_values(:, i_AR), 'LineWidth', 2, ...
%         'DisplayName', ['AR = ' num2str(AR_values(i_AR))], 'Color', colors(i_AR));
% end
% grid on;
% xlabel('AoA [deg]', 'Interpreter', 'latex', 'FontSize', 22);
% ylabel('$C_l$', 'Interpreter', 'latex', 'FontSize', 22);
% title('Lift Coefficient vs Angle of Attack');
% legend('Location', 'northwest', 'Interpreter', 'latex');

figure(1), clf;
hold on;
colors = ["#00a5cf", "#d1495b", "#3c1642", "#4f772d", "#edae49"];
xrange = [-7, 11];
yrange = [-0.5, 2];
fontsize_general = 18;


% enhance axes

xlim(xrange);
ylim(yrange);
plot(xlim, [0,0], 'LineWidth', 1, 'HandleVisibility', 'off', 'Color', 'k');
plot([0, 0], ylim, 'k', 'LineWidth', 1, 'HandleVisibility', 'off');
xticks(-6:2:10);

for i_AR = 1:length(AR_values)
    plot(alpha_data, Cl_values(:, i_AR), 'LineWidth', 2, ...
        'DisplayName', ['AR = ' num2str(AR_values(i_AR))], 'Color', colors(i_AR));
end
grid on

%title('Lift Coefficient vs Angle of Attack');
legend('Location', 'northwest', 'Interpreter', 'latex');
set(gca, 'Color', 'white', 'FontSize', fontsize_general ,'Xcolor', 'k', 'YColor', 'k', 'GridColor', [0.1, 0.1, 0.1], 'XMinorGrid', 'off', 'YMinorGrid', 'off');

xlabel('AoA [deg]', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$C_L$', 'Interpreter', 'latex', 'FontSize', 22);


%% Plot 2: CDi vs AoA
figure(2), clf;
hold on;
xrange = [-7, 11];
yrange = [-0.01, 0.1];

xticks(-6:2:10);
% enhance axes

xlim(xrange);
ylim(yrange);
plot(xlim, [0,0], 'LineWidth', 1, 'HandleVisibility', 'off', 'Color', 'k');
plot([0, 0], ylim, 'k', 'LineWidth', 1, 'HandleVisibility', 'off');


for i_AR = 1:length(AR_values)
    plot(alpha_data, CDi_values(:, i_AR), 'LineWidth', 2, ...
        'DisplayName', ['AR = ' num2str(AR_values(i_AR))], 'Color', colors(i_AR));
end

grid on;

%title('Induced Drag Coefficient vs Angle of Attack');
legend('Location', 'northwest', 'Interpreter','latex');


set(gca, 'Color', 'white', 'FontSize', fontsize_general ,'Xcolor', 'k', 'YColor', 'k', 'GridColor', [0.1, 0.1, 0.1], 'XMinorGrid', 'off', 'YMinorGrid', 'off');
xlabel('AoA [deg]', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$C_{D,i}$', 'Interpreter', 'latex', 'FontSize', 22);

%% Plot 3: ai vs AoA
figure(3), clf;
hold on;

xrange = [-7, 11];
yrange = [-1, 5];


% enhance axes

xlim(xrange);
xticks(-6:2:10);
ylim(yrange);
plot(xlim, [0,0], 'LineWidth', 1, 'HandleVisibility', 'off', 'Color', 'k');
plot([0, 0], ylim, 'k', 'LineWidth', 1, 'HandleVisibility', 'off');


for i_AR = 1:length(AR_values)
    plot(alpha_data, ai_values(:, i_AR), 'LineWidth', 2, ...
        'DisplayName', ['AR = ' num2str(AR_values(i_AR))], 'Color', colors(i_AR));
end
grid on;
legend('Location', 'northwest', 'Interpreter','latex');
set(gca, 'Color', 'white', 'FontSize', fontsize_general ,'Xcolor', 'k', 'YColor', 'k', 'GridColor', [0.1, 0.1, 0.1], 'XMinorGrid', 'off', 'YMinorGrid', 'off');
xlabel('AoA [deg]', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$\alpha_i$ [deg]', 'Interpreter', 'latex', 'FontSize', 22);


%% Plot 2: CDf vs AoA
figure(4), clf;
hold on;


xrange = [-7, 11];
yrange = [0.005, 0.012];


% enhance axes

xlim(xrange);
ylim(yrange);
xticks(-6:2:10);
yticks(0.005:0.001:0.012)

%plot(xlim, [0,0], 'LineWidth', 1, 'HandleVisibility', 'off', 'Color', 'k');
plot([0, 0], ylim, 'k', 'LineWidth', 1, 'HandleVisibility', 'off');

for i_AR = 1:length(AR_values)
    plot(alpha_data, CDf_values(:, i_AR), 'LineWidth', 2, ...
        'DisplayName', ['AR = ' num2str(AR_values(i_AR))], 'Color', colors(i_AR));
end
grid on;
legend('Location', 'northwest', 'Interpreter','latex');


set(gca, 'Color', 'white', 'FontSize', fontsize_general ,'Xcolor', 'k', 'YColor', 'k', 'GridColor', [0.1, 0.1, 0.1], 'XMinorGrid', 'off', 'YMinorGrid', 'off');
xlabel('AoA [deg]', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$C_D$', 'Interpreter', 'latex', 'FontSize', 22);