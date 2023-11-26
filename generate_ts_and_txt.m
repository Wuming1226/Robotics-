clear
clc

traject = hagongda([355 -138 486.96]', 90, 2);

figure(1)
plot3(traject(1,:), traject(2,:), traject(3,:), 'Linewidth', 2, 'Color', 'black', 'LineStyle', '-');
% hold on; 
grid on;

set(gca,'FontSize',24)
set(gcf,'unit','normalized','position', [0,0,0.465,0.8])
axis equal;

x_t = traject(1,:);
y_t = traject(2,:);
z_t = traject(3,:);
q_t = zeros(4, length(x_t));
q_t = q_t + [1 0 0 0]';

theta_start = [-20.9 14.5 127.5 180 -37.8 -150]'; 

theta_t = inverse_kinematics(x_t, y_t, z_t, q_t, theta_start);

theta1_ts = timeseries(theta_t(1, :), 0.001:0.001:length(theta_t)*0.001);
theta2_ts = timeseries(theta_t(2, :), 0.001:0.001:length(theta_t)*0.001);
theta3_ts = timeseries(theta_t(3, :), 0.001:0.001:length(theta_t)*0.001);
theta4_ts = timeseries(theta_t(4, :), 0.001:0.001:length(theta_t)*0.001);
theta5_ts = timeseries(theta_t(5, :), 0.001:0.001:length(theta_t)*0.001);
theta6_ts = timeseries(theta_t(6, :), 0.001:0.001:length(theta_t)*0.001);


% save data
fid = fopen('theta.txt', 'wt');
fprintf(fid, '%3.4f %3.4f %3.4f %3.4f %3.4f %3.4f\n', theta_t(1:6, :)*180/pi);
% fprintf(fid, '%3.4f %3.4f %3.4f %3.4f %3.4f %3.4f\n', theta_t(2:7, :)*180/pi);
fclose(fid);