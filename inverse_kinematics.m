% ABSTRACT: ������˶�ѧ������ؽڿռ�켣
%
% INPUT: x_t, y_t, z_t  ��ά�������У��ֱ�Ϊ 1��N ����
%        q_t            ��λ��Ԫ�����У�4��N ����
%        theta_start    ��ʼ�ؽ�λ�ã�������ѡ���˶�ѧ�Ľ�
%
% OUTPUT: theta_t       �ؽ�λ�����У�6��N ����
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