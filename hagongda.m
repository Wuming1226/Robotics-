% ABSTRACT: 规划书写汉字“哈工大”的三维坐标轨迹
%
% INPUT: position      	第一笔画的起点位置，3×1 向量
%        rotation       旋转角度
%		 scale          缩放倍数
%
% OUTPUT: traject	 	三维坐标轨迹点位序列，3×N 矩阵
%
function traject = hagongda(position, rotation, scale)

    max_depth = 5;
    bridge_height = 10;
    bridge_time = 1.5;

    % 笔画时间、高度
    % “哈”
    [x1, y1, z1] = generate_stroke('stroke\zhishu', [0 0] * scale, 10 * scale, max_depth);
    [tx1, ty1, tz1] = draw_stroke(x1, y1, z1);
    
    [x2, y2, z2] = generate_stroke('stroke\hengzhe', [0 0.75] * scale, 14 * scale, max_depth);
    [tx2, ty2, tz2] = draw_stroke(x2, y2, z2);
    
    [x3, y3, z3] = generate_stroke('stroke\zhongheng', [0.5 -10.5] * scale, 9 * scale, max_depth);
    [tx3, ty3, tz3] = draw_stroke(x3, y3, z3);
    
    [x4, y4, z4] = generate_stroke('stroke\changpie', [26 10] * scale, 40 * scale, max_depth);
    [tx4, ty4, tz4] = draw_stroke(x4, y4, z4);
    
    [x5, y5, z5] = generate_stroke('stroke\xiena', [23 5] * scale, 40 * scale, max_depth);
    [tx5, ty5, tz5] = draw_stroke(x5, y5, z5);
    
    [x6, y6, z6] = generate_stroke('stroke\zhongheng', [24 -15] * scale, 11 * scale, max_depth);
    [tx6, ty6, tz6] = draw_stroke(x6, y6, z6);
    
    [x7, y7, z7] = generate_stroke('stroke\zhishu', [22 -25] * scale, 14 * scale, max_depth);
    [tx7, ty7, tz7] = draw_stroke(x7, y7, z7);
    
    [x8, y8, z8] = generate_stroke('stroke\hengzhe', [22 -24.25] * scale, 22 * scale, max_depth);
    [tx8, ty8, tz8] = draw_stroke(x8, y8, z8);
    
    [x9, y9, z9] = generate_stroke('stroke\zhongheng', [22.5 -36.5] * scale, 15 * scale, max_depth);
    [tx9, ty9, tz9] = draw_stroke(x9, y9, z9);
    
    % “工”
    [x10, y10, z10] = generate_stroke('stroke\changheng', [77 00] * scale, 25 * scale, max_depth);
    [tx10, ty10, tz10] = draw_stroke(x10, y10, z10);
    
    [x11, y11, z11] = generate_stroke('stroke\zhishu', [89.5 -2] * scale, 25 * scale, max_depth);
    [tx11, ty11, tz11] = draw_stroke(x11, y11, z11);
    
    [x12, y12, z12] = generate_stroke('stroke\changheng', [67 -30] * scale, 45 * scale, max_depth);
    [tx12, ty12, tz12] = draw_stroke(x12, y12, z12);
    
    % “大”
    [x13, y13, z13] = generate_stroke('stroke\changheng', [135 -5] * scale, 28 * scale, max_depth);
    [tx13, ty13, tz13] = draw_stroke(x13, y13, z13);
    
    [x14, y14, z14] = generate_stroke('stroke\hupie', [147.5 10] * scale, 50 * scale, max_depth);
    [tx14, ty14, tz14] = draw_stroke(x14, y14, z14);
    
    [x15, y15, z15] = generate_stroke('stroke\xiena', [147.5 -12] * scale, 40 * scale, max_depth);
    [tx15, ty15, tz15] = draw_stroke(x15, y15, z15);
    
    % 连接段
    [bridge_x1, bridge_y1, bridge_z1] = bridge_connect([tx1(end) ty1(end) tz1(end)]', [tx2(1) ty2(1) tz2(1)]', bridge_height, bridge_time);
    [bridge_x2, bridge_y2, bridge_z2] = bridge_connect([tx2(end) ty2(end) tz2(end)]', [tx3(1) ty3(1) tz3(1)]', bridge_height, bridge_time);
    [bridge_x3, bridge_y3, bridge_z3] = bridge_connect([tx3(end) ty3(end) tz3(end)]', [tx4(1) ty4(1) tz4(1)]', bridge_height, bridge_time);
    [bridge_x4, bridge_y4, bridge_z4] = bridge_connect([tx4(end) ty4(end) tz4(end)]', [tx5(1) ty5(1) tz5(1)]', bridge_height, bridge_time);
    [bridge_x5, bridge_y5, bridge_z5] = bridge_connect([tx5(end) ty5(end) tz5(end)]', [tx6(1) ty6(1) tz6(1)]', bridge_height, bridge_time);
    [bridge_x6, bridge_y6, bridge_z6] = bridge_connect([tx6(end) ty6(end) tz6(end)]', [tx7(1) ty7(1) tz7(1)]', bridge_height, bridge_time);
    [bridge_x7, bridge_y7, bridge_z7] = bridge_connect([tx7(end) ty7(end) tz7(end)]', [tx8(1) ty8(1) tz8(1)]', bridge_height, bridge_time);
    [bridge_x8, bridge_y8, bridge_z8] = bridge_connect([tx8(end) ty8(end) tz8(end)]', [tx9(1) ty9(1) tz9(1)]', bridge_height, bridge_time);
    [bridge_x9, bridge_y9, bridge_z9] = bridge_connect([tx9(end) ty9(end) tz9(end)]', [tx10(1) ty10(1) tz10(1)]', bridge_height, bridge_time);
    [bridge_x10, bridge_y10, bridge_z10] = bridge_connect([tx10(end) ty10(end) tz10(end)]', [tx11(1) ty11(1) tz11(1)]', bridge_height, bridge_time);
    [bridge_x11, bridge_y11, bridge_z11] = bridge_connect([tx11(end) ty11(end) tz11(end)]', [tx12(1) ty12(1) tz12(1)]', bridge_height, bridge_time);
    [bridge_x12, bridge_y12, bridge_z12] = bridge_connect([tx12(end) ty12(end) tz12(end)]', [tx13(1) ty13(1) tz13(1)]', bridge_height, bridge_time);
    [bridge_x13, bridge_y13, bridge_z13] = bridge_connect([tx13(end) ty13(end) tz13(end)]', [tx14(1) ty14(1) tz14(1)]', bridge_height, bridge_time);
    [bridge_x14, bridge_y14, bridge_z14] = bridge_connect([tx14(end) ty14(end) tz14(end)]', [tx15(1) ty15(1) tz15(1)]', bridge_height, bridge_time);

    % 拼接
    traject = [tx1 bridge_x1 tx2 bridge_x2 tx3 bridge_x3 tx4 bridge_x4 tx5 bridge_x5 tx6 bridge_x6 tx7 bridge_x7 tx8 bridge_x8 tx9 bridge_x9 tx10 bridge_x10 tx11 bridge_x11 tx12 bridge_x12 tx13 bridge_x13 tx14 bridge_x14 tx15;
        ty1 bridge_y1 ty2 bridge_y2 ty3 bridge_y3 ty4 bridge_y4 ty5 bridge_y5 ty6 bridge_y6 ty7 bridge_y7 ty8 bridge_y8 ty9 bridge_y9 ty10 bridge_y10 ty11 bridge_y11 ty12 bridge_y12 ty13 bridge_y13 ty14 bridge_y14 ty15;
        tz1 bridge_z1 tz2 bridge_z2 tz3 bridge_z3 tz4 bridge_z4 tz5 bridge_z5 tz6 bridge_z6 tz7 bridge_z7 tz8 bridge_z8 tz9 bridge_z9 tz10 bridge_z10 tz11 bridge_z11 tz12 bridge_z12 tz13 bridge_z13 tz14 bridge_z14 tz15];
    
    % 旋转
    R = [cos(rotation*pi/180) -sin(rotation*pi/180) 0; sin(rotation*pi/180) cos(rotation*pi/180) 0; 0 0 1];
    traject = R * traject;
    
    % 位置偏置
    traject = traject + position;
end


% ABSTRACT: 规划笔画的三维坐标轨迹
%
% INPUT: x, y, z      	笔画路径点序列，分别为 1×N 矩阵
%
% OUTPUT: tx, ty, tz    三维坐标轨迹点位序列，分别为 1×N 矩阵
%
function [tx, ty, tz] = draw_stroke(x, y, z)
    ddtheta = 1200;
    interval = 0.1;
    step = 0.001;
    [tx, ~, ~] = LFPB(x, interval * ones(1, length(x)-1), ddtheta, step, 0);
    [ty, ~, ~] = LFPB(y, interval * ones(1, length(y)-1), ddtheta, step, 0);
    [tz, ~, ~] = LFPB(z, interval * ones(1, length(z)-1), ddtheta, step, 0);
end


% ABSTRACT: 规划笔画间连接段三维坐标轨迹
%
% INPUT: start          连接段起点，上一笔画的终点，3×1 向量
%		 term           连接段终点，下一笔画的起点，3×1 向量
%		 height         连接段高度
%		 time           连接段时间
%
% OUTPUT: bridge_x, bridge_y, bridge_z	 	三维坐标轨迹点位序列，分别为 1×N 矩阵
%
function [bridge_x, bridge_y, bridge_z] = bridge_connect(start, term, height, time)
    waypoint1 = start + [0 0 height]';
    waypoint2 = term + [0 0 height]';
    
    ddtheta = 300;
    step = 0.001;
    [bridge_x, ~, ~] = LFPB([start(1) waypoint1(1) waypoint2(1) term(1)], [time/3 time/3 time/3], ddtheta, step, 0);
    [bridge_y, ~, ~] = LFPB([start(2) waypoint1(2) waypoint2(2) term(2)], [time/3 time/3 time/3], ddtheta, step, 0);
    [bridge_z, ~, ~] = LFPB([start(3) waypoint1(3) waypoint2(3) term(3)], [time/3 time/3 time/3], ddtheta, step, 0); 
end