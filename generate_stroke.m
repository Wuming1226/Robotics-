% ABSTRACT: ���ɱʻ���ά��������
%
% INPUT: stroke         �ʻ����ƣ��ַ���
%        position      	��һ�ʻ���ƽ������λ�ã�2��1 ����
%        size           �ʻ��ߴ�
%		 max_depth      �����ѹ���
%
% OUTPUT: x_list, y_list, z_list    ��ά����켣��λ���У��ֱ�Ϊ 1��N ����
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