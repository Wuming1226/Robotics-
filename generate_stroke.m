% ABSTRACT: 生成笔画三维坐标序列
%
% INPUT: stroke         笔画名称，字符串
%        position      	第一笔画在平面的起点位置，2×1 向量
%        size           笔画尺寸
%		 max_depth      最大下压深度
%
% OUTPUT: x_list, y_list, z_list    三维坐标轨迹点位序列，分别为 1×N 矩阵
%
function [x_list, y_list, z_list] = generate_stroke(stroke, position, size, max_depth)
    x_list = load([stroke, '_x.txt']);
    y_list = load([stroke, '_y.txt']);
    
    len = sqrt(([x_list(1); y_list(1)] - [x_list(end); y_list(end)])'* ([x_list(1); y_list(1)] - [x_list(end); y_list(end)]));
    x_list = x_list .* (size / len);
    y_list = y_list .* (size / len);

    x_list = x_list + position(1);
    y_list = y_list + position(2);
    
    N = length(x_list);
    z_list = -max_depth * ones(1, N);
    for i = 1:3
        z_list(i) = -max_depth / 4 * i;
    end
    for i = N-2:N
        z_list(i) = -max_depth / 4 * (N - i + 1);
    end
end