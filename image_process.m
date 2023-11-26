clear
clc

src = imread('stroke_image/changheng.png');

thresh = graythresh(src);
bin = im2bw(src, thresh);

bin = imcomplement(bin);

figure(1)
imshow(bin)

% mask = [0 0 0 0 0; 0 0 1 0 0; 0 1 1 0 0; 0 0 1 0 0; 0 0 0 0 0];
% erd = imerode(bin, mask);
% mask = [0 0 0 0 0; 0 0 0 0 0; 0 1 1 1 0; 0 0 0 0 0; 0 0 0 0 0];
% dia = imdilate(erd, mask);
% 
mask = [0 0 0 0 0; 0 1 1 1 0; 0 1 1 1 0; 0 1 1 1 0; 0 0 0 0 0];
% mask = [1 1 1 1 1; 1 1 1 1 1; 1 1 1 1 1; 1 1 1 1 1; 1 1 1 1 1];
dia = imdilate(bin, mask);
mask = [0 0 0 0 0; 0 0 1 0 0; 0 1 1 1 0; 0 0 1 0 0; 0 0 0 0 0];
% mask = [0 0 0 0 0; 0 0 0 0 0; 0 0 1 0 0; 0 0 0 0 0; 0 0 0 0 0];
erd = imerode(dia, mask);

figure(2)
imshow(dia)
figure(3)
imshow(erd)

% L = bwlabel(dia, 8);
L = bwlabel(erd, 8);

[m, n] = size(L);
left = zeros(m, n);
for r = 1:m
    for c = 1:n
        if L(r, c) == 1
            left(r, c) = 1;
        end
    end
end


start_r = 83;
start_c = 61;
dst = zeros(m, n);
dst(start_r, start_c) = 1;
ldst = dst;
left(start_r, start_c) = 0;
s = 1;
while sum(sum(left)) ~= 0
    s = s + 1;
    for r = 2:m-1
        for c = 2:n-1
            if left(r, c) == 1
                if ldst(r-1, c) ~= 0 || ldst(r+1, c) ~= 0 || ldst(r, c-1) ~= 0 || ldst(r, c+1) ~= 0 || ldst(r-1, c-1) ~= 0 || ldst(r+1, c-1) ~= 0 || ldst(r-1, c-1) ~= 0 || ldst(r+1, c+1) ~= 0
                    dst(r, c) = s;
                    left(r, c) = 0;
                end
                
            end
        end
    end
    ldst = dst;
end

step = 10;
x_list = zeros(1, floor(max(max(dst)) / step)+1);
y_list = zeros(1, floor(max(max(dst)) / step)+1);
x_list(1) = start_c;
y_list(1) = n - start_r;
for i=1:1:(floor(max(max(dst)) / step))
    cx = 0;
    cy = 0;
    cnt = 0;
    for r = 2:m-1
        for c = 2:n-1
            if dst(r, c) == i * step
                cx = cx + c;
                cy = cy + n - r;
                cnt = cnt + 1;
            end
        end
    end 
    x_list(i+1) = cx / cnt;
    y_list(i+1) = cy / cnt;
end

x_list = x_list - x_list(1);
y_list = y_list - y_list(1);

figure(1)
set(gcf,'unit','normalized','position', [0,0,0.465,0.8])

plot(x_list, y_list,  'Linewidth', 2, 'Color', 'black', 'LineStyle', '-'); 
grid on;
axis equal;

% save data
% fid = fopen('stroke/changheng_x.txt', 'wt');
% fprintf(fid, '%3.4f\n', x_list);
% fid = fopen('stroke/changheng_y.txt', 'wt');
% fprintf(fid, '%3.4f\n', y_list);
% fclose(fid);

