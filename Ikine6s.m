% ABSTRACT: 求解六轴机器人逆运动学所有解
%
% INPUT: T         目标位姿，刚体变换矩阵
%        config    六轴机器人尺寸信息
%
% OUTPUT: theta    逆运动学解，6×N 矩阵
%
function theta = Ikine6s(T, config)

    theta = [];
    
    Xi = [0 -config(1) -config(1)-config(2) 0 -config(1)-config(2)-config(3) 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0;
    0 1 1 0 1 0;
    1 0 0 1 0 1];

    p1 = [0 0 config(1)+config(2)+config(3)]';
    p2 = [0 0 config(1)]';
    p3 = [0 0 config(1)]';
    p4 = [config(1) 0 0]';

    g0 = [-1 0 0 0; 0 -1 0 0; 0 0 1 config(1)+config(2)+config(3)+config(4); 0 0 0 1];

    g1 = T / g0;
    g1p1 = g1*[p1;1];
    g1p1 = g1p1(1:3);
    theta3 = subprob3(p1, p2, sqrt((g1p1 - p2)'* (g1p1 - p2)), Xi(4:6, 3), [0 0 config(1)+config(2)]');
   
    for index1=1:length(theta3)
        e3p1 = exp_xi(Xi(:, 3), theta3(index1)) * [p1;1];
        e3p1 = e3p1(1:3);
        [theta1, theta2] = subprob2(e3p1, g1p1, Xi(4:6, 1), Xi(4:6, 2), p2);
        
        for index2=1:length(theta1)
            g2 = (exp_xi(Xi(:, 1), theta1(index2)) * exp_xi(Xi(:, 2), theta2(index2)) * exp_xi(Xi(:, 3), theta3(index1))) \ g1;
            g2p3 = g2 * [p3;1];
            g2p3 = g2p3(1:3);
            [theta4, theta5] = subprob2(p3, g2p3, Xi(4:6, 4), Xi(4:6, 5), p1);
            
            for index3=1:length(theta4)
                g3 = (exp_xi(Xi(:, 4), theta4(index3)) * exp_xi(Xi(:, 5), theta5(index3))) \ g2;
                g3p4 = g3 * [p4; 1];
                g3p4 = g3p4(1:3);
                theta6 = subprob1(p4, g3p4, Xi(4:6, 6), p1);
                
                for index4=1:length(theta6)
                    theta = [theta, [theta1(index2); theta2(index2); theta3(index1); theta4(index3); theta5(index3); theta6(index4)]];
                end
            end
        end
    end
    
    for r = 1:6
        for c = 1:length(theta)
            if theta(r, c) <= -pi + 0.01
                theta(r, c) = theta(r, c) + 2*pi;
            end
        end
    end
    
end
   

function wh = hat(w)
    wh = [0 -w(3) w(2); w(3) 0 -w(1); -w(2) w(1) 0];
end


function R = exp_w(w, theta)
    R = eye(3) + hat(w) * sin(theta) + hat(w)^2 * (1 - cos(theta));
end


function g = exp_xi(Xi, theta)
    v = Xi(1:3);
    w = Xi(4:6);
    R = exp_w(w, theta);
    p = (eye(3) - R) * hat(w) * v + w * w' * v * theta;
    
    g = [R p; 0 0 0 1];
end


function theta = subprob1(p, q, w, r)
    u = p - r;
    v = q - r;
    u_proj = (eye(3) - w * w') * u;
    v_proj = (eye(3) - w * w') * v;
    if (w'* u) - (w'* v) > 0.0001 || (w'* u) - (w'* v) < -0.0001 || (u_proj'* u_proj) - (v_proj'* v_proj) > 0.0001 || (u_proj'* u_proj) - (v_proj'* v_proj) < -0.0001
        theta = 999;
    else
        theta = atan2(w'* cross(u_proj, v_proj), u_proj'* v_proj);
    end
end


function [theta1, theta2] = subprob2(p, q, w1, w2, r)
    if cross(w1, w2) == 0
        theta1 = subprob1(p, q, w1, r);
        theta2 = 0;
    else
        u = p - r;
        v = q - r;
        alpha = ((w1'* w2) * w2'* u - w1'* v) / ((w1'* w2)' * (w1'* w2) - 1);
        beta = ((w1'* w2) * w1'* v - w2'* u) / ((w1'* w2)'* (w1'* w2) - 1);
        gamma_sq = (u'* u - alpha^2 - beta^2 - 2 * alpha * beta * w1'* w2) / (cross(w1, w2)'* cross(w1, w2));
        if gamma_sq < 0
            error('No solution for Subproblem 2.')
        elseif gamma_sq == 0
            z = alpha * w1 + beta * w2;
            c = r + z;
            theta1 = subprob1(c, q, w1, r);
            theta2 = subprob1(p, c, w2, r);
        else
           gamma1 = sqrt(gamma_sq);
           gamma2 = -sqrt(gamma_sq);
           z1 = alpha * w1 + beta * w2 + gamma1 * cross(w1, w2);
           z2 = alpha * w1 + beta * w2 + gamma2 * cross(w1, w2);
           c1 = r + z1;
           c2 = r + z2;
           theta11 = subprob1(c1, q, w1, r);
           theta21 = subprob1(p, c1, w2, r);
           theta12 = subprob1(c2, q, w1, r);
           theta22 = subprob1(p, c2, w2, r);
           theta1 = [theta11 theta12];
           theta2 = [theta21 theta22];
        end
    end
end


function theta = subprob3(p, q, delta, w, r)
    u = p - r;
    v = q - r;
    u_proj = u - w * w'* u;
    v_proj = v - w * w'* v;
    if delta^2 < (w'*(p - q))'* (w'*(p - q))
        error('No solution for Subproblem 3')
    else
        delta_proj_sq = delta^2 - (w'*(p - q))'* (w'*(p - q));
        cs = (u_proj'* u_proj + v_proj'* v_proj - delta_proj_sq) / (2 * sqrt(u_proj'* u_proj) * sqrt(v_proj'* v_proj));
        if cs > 1 || cs < -1
            error('No solution for Subproblem 3.')
        elseif cs == 1
            theta = atan2(w'* cross(u_proj, v_proj), u_proj'* v_proj);
        elseif cs == -1
            theta = atan2(w'* cross(u_proj, v_proj), u_proj'* v_proj) - pi;
        else
            theta1 = atan2(w'* cross(u_proj, v_proj), u_proj'* v_proj) - acos(cs);
            theta2 = atan2(w'* cross(u_proj, v_proj), u_proj'* v_proj) + acos(cs) - 2*pi;
            theta = [theta1, theta2];
        end
    end 
end