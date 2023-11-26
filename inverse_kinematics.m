% ABSTRACT: 求解逆运动学、计算关节空间轨迹
%
% INPUT: x_t, y_t, z_t  三维坐标序列，分别为 1×N 矩阵
%        q_t            单位四元数序列，4×N 矩阵
%        theta_start    初始关节位置，用于挑选逆运动学的解
%
% OUTPUT: theta_t       关节位置序列，6×N 矩阵
%
function theta_t = inverse_kinematics(x_t, y_t, z_t, q_t, theta_start)

    theta0 = theta_start;
    theta_t = zeros(6, length(x_t));
    
    for i=1:length(x_t)
        p = [x_t(i) y_t(i) z_t(i)]';
        R = eye(3) + 2 * q_t(1, i) * hat(q_t(2:4, i)) + 2 * hat(q_t(2:4, i))^2;
        
        R0 = exp_w([1 0 0], pi);
        config = [491 450 450 84]';
        
        T = [R0 * R p;
            0 0 0 1];
        
        theta_r = Ikine6s(T, config);
        [~, index] = min(sum(abs(theta_r - theta0), 1));
        theta_r = theta_r(1:6, index);

        theta_t(:, i) = theta_r;
        theta0 = theta_r;
    end
end

function wh = hat(w)
    wh = [0 -w(3) w(2); w(3) 0 -w(1); -w(2) w(1) 0];
end

function R = exp_w(w, theta)
    R = eye(3) + hat(w) * sin(theta) + hat(w)^2 * (1 - cos(theta));
end