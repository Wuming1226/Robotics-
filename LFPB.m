% ABSTRACT: 使用 LFPB 方法，根据给定路径点和对应转移时间，生成平滑位置轨迹
%
% INPUT: theta      	路径点序列，1×N 矩阵
%		 t_d            转移时间序列，1×N-1 矩阵
%        ddtheta_m      最大加速度的绝对值
%        step           时间步长
%        t_0            初始时间
%
% OUTPUT: theta_t	 	位置序列，1×N 矩阵
%         dtheta_t      时间-速度序列，2×N 矩阵
%         ddtheta_t     时间-加速度序列，2×N 矩阵
%
function [theta_t, dtheta_t, ddtheta_t] = LFPB(theta, t_d, ddtheta_m, step, t0)
    
    N = length(theta);      % 路径点序列长度
    theta_t = zeros(2, round(sum(t_d) / step));
    dtheta_t = zeros(2, round(sum(t_d) / step));
    ddtheta_t = zeros(2, round(sum(t_d) / step));
    
    dtheta = (theta(2:N) - theta(1:N-1)) ./ t_d;    % 直线轨迹速度序列
    
    t_b = zeros(1, N);      % 抛物线轨迹时间序列
    t_l = zeros(1, N-1);    % 直线轨迹时间序列
    
    % 初始段参数
    if abs(theta(2) - theta(1)) < 0.00001
        ddtheta_1 = 0;
        t_b(1) = 0;
    else
        ddtheta_1 = sign(theta(2) - theta(1)) * ddtheta_m;
        t_b(1) = t_d(1) - sqrt(t_d(1)^2 - 2 * (theta(2) - theta(1)) / ddtheta_1);
    end
    dtheta(1) = (theta(2) - theta(1)) / (t_d(1) - t_b(1) / 2);
    
    % 结束段参数
    if abs(theta(N) - theta(N-1)) < 0.00001
        ddtheta_N = 0;
        t_b(N) = 0;
    else
        ddtheta_N = sign(theta(N) - theta(N-1)) * ddtheta_m;
        t_b(N) = t_d(N-1) - sqrt(t_d(N-1)^2 - 2 * (theta(N) - theta(N-1)) / ddtheta_N);
    end
    dtheta(N-1) = (theta(N) - theta(N-1)) / (t_d(N-1) - t_b(N) / 2);
    
    % 生成初始段抛物线轨迹
    for i=1:round((t_b(1))/step)
        t = i * step;
        theta_t(1, i) = t;
        dtheta_t(1, i) = t;
        ddtheta_t(1, i) = t;
        theta_t(2, i) = theta(1) + 0.5 * ddtheta_1 * t^2;
        dtheta_t(2, i) = ddtheta_1 * t;
        ddtheta_t(2, i) = ddtheta_1;
    end
    linear_start = theta(1) + 0.5 * ddtheta_1 * t_b(1)^2;   % 直线轨迹起始位置
    
    % 生成中间段直线-抛物线轨迹组
    for seg = 1:N-1
        t_start = sum(t_b(1:seg)) + sum(t_l(1:seg));    % 直线轨迹初始时间
        
        if seg < N-1    % 结束段不更新 ddtheta_b 和 t_b
            if abs(dtheta(seg+1) - dtheta(seg)) < 0.00001
                ddtheta_b = 0;
                t_b(seg+1) = 0;
            else
                ddtheta_b = sign(dtheta(seg+1) - dtheta(seg)) * ddtheta_m;
                t_b(seg+1) = (dtheta(seg+1) - dtheta(seg)) / ddtheta_b;
            end
        end
        if seg == 1     % 初始段计算公式特殊
            t_l(1) = t_d(1) - t_b(1) - t_b(2) / 2;
        else
            t_l(seg) = t_d(seg) - t_b(seg) / 2 - t_b(seg+1) / 2;
        end
        
        % 生成直线轨迹
        for i=round(t_start/step)+1: round((t_start + t_l(seg))/step)
            t = i * step;
            theta_t(1, i) = t;
            dtheta_t(1, i) = t;
            ddtheta_t(1, i) = t;
            theta_t(2, i) = linear_start + dtheta(seg) * (t - t_start);
            dtheta_t(2, i) = dtheta(seg);
            ddtheta_t(2, i) = 0;
        end
        blend_start = linear_start + dtheta(seg) * t_l(seg);    % 抛物线轨迹起始位置
        if seg < N-1 && t_b(seg+1) >= step     % 若抛物线轨迹时间大于step，生成抛物线轨迹（倒数第二段只生成直线轨迹）
            for i=round((t_start + t_l(seg))/step)+1: round((t_start + t_l(seg) + t_b(seg+1))/step)
                t = i * step;
                theta_t(1, i) = t;
                dtheta_t(1, i) = t;
                ddtheta_t(1, i) = t;
                % 抛物线轨迹参数计算
                m = t_start + t_l(seg) - dtheta(seg) / ddtheta_b;
                b = blend_start - 0.5 * dtheta(seg)^2 / ddtheta_b;
                theta_t(2, i) = b + 0.5 * ddtheta_b * (t - m)^2;
                dtheta_t(2, i) = ddtheta_b * (t - m);
                ddtheta_t(2, i) = ddtheta_b;
            end
            linear_start = b + 0.5 * ddtheta_b * (t_start + t_l(seg) + t_b(seg+1) - m)^2;   % 直线轨迹起始位置
        else
            linear_start = blend_start;
        end 
        
    end

    % 生成结束段抛物线
    for i=round((sum(t_d)-t_b(N))/step)+1: round(sum(t_d)/step)
        t = i * step;
        theta_t(1, i) = t;
        dtheta_t(1, i) = t;
        ddtheta_t(1, i) = t;
        theta_t(2, i) = theta(N) - 0.5 * ddtheta_N  * (sum(t_d) - t)^2;
        dtheta_t(2, i) = ddtheta_N * (sum(t_d) - t);
        ddtheta_t(2, i) = -ddtheta_N;
    end
    
    % 加入初始位置
    theta_t = [[0; theta(1)] theta_t]; 
    dtheta_t = [[0; 0] dtheta_t]; 
    ddtheta_t = [[0; 0] ddtheta_t]; 
    
    % 初始时间偏置
    theta_t(1, 1:round(sum(t_d)/step)) = theta_t(1, 1:round(sum(t_d)/step)) + t0;
    dtheta_t(1, 1:round(sum(t_d)/step)) = dtheta_t(1, 1:round(sum(t_d)/step)) + t0;
    ddtheta_t(1, 1:round(sum(t_d)/step)) = ddtheta_t(1, 1:round(sum(t_d)/step)) + t0;
    
    theta_t(1, :) = [];
end
