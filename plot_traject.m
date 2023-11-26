clear
clc

traject = hagongda([355 -138 486.96]', 90, 2);

% traject = traject(:, 1:206014);

figure(1)
plot3(traject(1,:), traject(2,:), traject(3,:), 'Linewidth', 2, 'Color', 'black', 'LineStyle', '-');
% plot(traject(1,:), 'Linewidth', 2, 'Color', 'black', 'LineStyle', '-');
% hold on; 
grid on;

set(gca,'FontSize',24)
set(gcf,'unit','normalized','position', [0,0,0.465,0.8])
axis equal;

% x_t = traject(1,:);
% y_t = traject(2,:);
% z_t = traject(3,:);
% q_t = zeros(4, length(x_t));
% q_t = q_t + [1 0 0 0]';
% 
% theta_start = [-20.9 14.5 127.5 180 -37.8 -150]'; 
% 
% theta_t = inverse_kinematics(x_t, y_t, z_t, q_t, theta_start);