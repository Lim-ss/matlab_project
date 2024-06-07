% 《应力波基础讲义》
% 弹塑性加载波在固壁端的反射

template = struct('x', 0, 'y', 0);

% 绘制v-sigma图
subplot(1, 2, 1);
hold on

% 绘制弹性波AB
plot([0, 5], [0, 5], "black", 'LineWidth', 1);

% 绘制区域BCED
matrix_1 = repmat(template, 6, 16);
matrix_1(1, 1).x = 5;
matrix_1(1, 1).y = 5;

for j = 2:16
    matrix_1(1, j).x = matrix_1(1, j - 1).x + 1;
    matrix_1(1, j).y = matrix_1(1, j - 1).y + slope(matrix_1(1, j - 1).y);
    plot([matrix_1(1, j - 1).x, matrix_1(1, j).x], [matrix_1(1, j - 1).y, matrix_1(1, j).y], "yellow", 'LineWidth', 1);
end

for i = 2:6

    matrix_1(i, 1).x = matrix_1(i - 1, 1).x - 1;
    matrix_1(i, 1).y = matrix_1(i - 1, 1).y + slope(matrix_1(i - 1, 1).y);
    plot([matrix_1(i - 1, 1).x, matrix_1(i, 1).x], [matrix_1(i - 1, 1).y, matrix_1(i, 1).y], "yellow", 'LineWidth', 1);

    for j = 2:16
        matrix_1(i, j).x = matrix_1(i, j - 1).x + 1;
        matrix_1(i, j).y = matrix_1(i, j - 1).y + slope(matrix_1(i, j - 1).y);
        plot([matrix_1(i, j - 1).x, matrix_1(i, j).x], [matrix_1(i, j - 1).y, matrix_1(i, j).y], "yellow", 'LineWidth', 1);
        plot([matrix_1(i - 1, j).x, matrix_1(i, j).x], [matrix_1(i - 1, j).y, matrix_1(i, j).y], "yellow", 'LineWidth', 1);
    end

end

% 绘制区域CEG
matrix_2 = repmat(template, 6, 6);

for i = 1:6
    matrix_2(i, 1).x = matrix_1(7 - i, 16).x;
    matrix_2(i, 1).y = matrix_1(7 - i, 16).y;
end

for j = 2:6
    matrix_2(1, j).x = matrix_2(1, j - 1).x + 1;
    matrix_2(1, j).y = matrix_2(1, j - 1).y + slope(matrix_2(1, j - 1).y);
    plot([matrix_2(1, j).x, matrix_2(1, j - 1).x], [matrix_2(1, j).y, matrix_2(1, j - 1).y], "red", 'LineWidth', 1);
end

for i = 2:6

    for j = 2:7 - i
        matrix_2(i, j).x = matrix_2(i, j - 1).x + 1;
        matrix_2(i, j).y = matrix_2(i, j - 1).y + slope(matrix_2(i, j - 1).y);
        plot([matrix_2(i, j - 1).x, matrix_2(i, j).x], [matrix_2(i, j - 1).y, matrix_2(i, j).y], "red", 'LineWidth', 1);
        plot([matrix_2(i - 1, j).x, matrix_2(i, j).x], [matrix_2(i - 1, j).y, matrix_2(i, j).y], "red", 'LineWidth', 1);
    end

end

% 绘制区域DEF
matrix_3 = repmat(template, 16, 16);

for j = 1:16
    matrix_3(1, j).x = matrix_1(6, 17 - j).x;
    matrix_3(1, j).y = matrix_1(6, 17 - j).y;
end

for i = 2:16
    matrix_3(i, 1).x = matrix_3(i - 1, 1).x - 1;
    matrix_3(i, 1).y = matrix_3(i - 1, 1).y + slope(matrix_3(i - 1, 1).y);
    plot([matrix_3(i, 1).x, matrix_3(i - 1, 1).x], [matrix_3(i, 1).y, matrix_3(i - 1, 1).y], "green", 'LineWidth', 1);
end

for j = 2:16

    for i = 2:17 - j
        matrix_3(i, j).x = matrix_3(i - 1, j).x - 1;
        matrix_3(i, j).y = matrix_3(i - 1, j).y + slope(matrix_3(i-1, j).y);
        plot([matrix_3(i, j).x, matrix_3(i - 1, j).x], [matrix_3(i, j).y, matrix_3(i - 1, j).y], "green", 'LineWidth', 1);
        plot([matrix_3(i, j).x, matrix_3(i, j - 1).x], [matrix_3(i, j).y, matrix_3(i, j - 1).y], "green", 'LineWidth', 1);
    end

end

% 绘制区域EFHG
matrix_4 = repmat(template, 16, 6);

for i = 1:16
    matrix_4(i, 1).x = matrix_3(i, 1).x;
    matrix_4(i, 1).y = matrix_3(i, 1).y;
end

for j = 1:6
    matrix_4(1, j).x = matrix_2(1, j).x;
    matrix_4(1, j).y = matrix_2(1, j).y;
end

for i = 2:16

    for j = 2:6
        matrix_4(i, j).x = matrix_4(i, j-1).x + 1;
        matrix_4(i, j).y = matrix_4(i, j-1).y + slope(matrix_4(i, j-1).y);
        plot([matrix_4(i, j).x, matrix_4(i - 1, j).x], [matrix_4(i, j).y, matrix_4(i - 1, j).y], "blue", 'LineWidth', 1);
        plot([matrix_4(i, j).x, matrix_4(i, j - 1).x], [matrix_4(i, j).y, matrix_4(i, j - 1).y], "blue", 'LineWidth', 1);
    end

end

% 绘制坐标轴和标记
plot([0, 0], [0, 25], "black", 'LineWidth', 1);
plot([0, 21], [0, 0], "black", 'LineWidth', 1);
plot([20, 20], [0, 25], "black--", 'LineWidth', 1);
plot([0, 5], [5, 5], "black--", 'LineWidth', 1);

text(1, 4, "\sigma_Y", 'FontSize', 12, 'Color', 'black', 'HorizontalAlignment', 'center');
text(1, 1, "A", 'FontSize', 12, 'Color', 'black', 'HorizontalAlignment', 'center');
text(6, 5, "B", 'FontSize', 12, 'Color', 'black', 'HorizontalAlignment', 'center');
text(1, matrix_1(6,1).y-1, "D", 'FontSize', 12, 'Color', 'black', 'HorizontalAlignment', 'center');
text(19, matrix_1(1,16).y, "C", 'FontSize', 12, 'Color', 'black', 'HorizontalAlignment', 'center');
text(matrix_1(6,16).x, matrix_1(6,16).y-1, "E", 'FontSize', 12, 'Color', 'black', 'HorizontalAlignment', 'center');
text(matrix_2(1,6).x-1, matrix_2(1,6).y+1, "G", 'FontSize', 12, 'Color', 'black', 'HorizontalAlignment', 'center');
text(matrix_3(16,1).x+1, matrix_3(16,1).y-1, "F", 'FontSize', 12, 'Color', 'black', 'HorizontalAlignment', 'center');
text(matrix_4(16,6).x, matrix_4(16,6).y+1, "H", 'FontSize', 12, 'Color', 'black', 'HorizontalAlignment', 'center');
text(21, 1, "v*", 'FontSize', 12, 'Color', 'black', 'HorizontalAlignment', 'center');

% 绘制X-t图
subplot(1, 2, 2);
hold on

plot([0, 0], [0, 20], "black", 'LineWidth', 1);
plot([0, 11], [0, 0], "black", 'LineWidth', 1);
plot([10, 10], [0, 20], "black--", 'LineWidth', 1);

plot([10, 0], [0, 5], "black", 'LineWidth', 1);

plot([0, 10], [5, 10], "blue", 'LineWidth', 1);
plot([0, 10], [5, 10.5], "blue", 'LineWidth', 1);
plot([0, 10], [5, 11], "blue", 'LineWidth', 1);
plot([0, 10], [5, 11.5], "blue", 'LineWidth', 1);
plot([10, 2], [10, 14], "blue", 'LineWidth', 1);
plot([10, 2], [10.5, 15], "blue", 'LineWidth', 1);
plot([10, 2], [11, 16], "blue", 'LineWidth', 1);
plot([10, 2], [11.5, 17], "blue", 'LineWidth', 1);

plot([10, 0], [0, 8], "red", 'LineWidth', 1);
plot([10, 0], [0, 8.5], "red", 'LineWidth', 1);
plot([10, 0], [0, 9], "red", 'LineWidth', 1);
plot([10, 0], [0, 9.5], "red", 'LineWidth', 1);
plot([0, 8], [8, 14], "red", 'LineWidth', 1);
plot([0, 8], [8.5, 15], "red", 'LineWidth', 1);
plot([0, 8], [9, 16], "red", 'LineWidth', 1);
plot([0, 8], [9.5, 17], "red", 'LineWidth', 1);

text(3, 2, "A", 'FontSize', 12, 'Color', 'black', 'HorizontalAlignment', 'center');
text(3, 4.5, "B", 'FontSize', 12, 'Color', 'black', 'HorizontalAlignment', 'center');
text(1, 6.5, "D", 'FontSize', 12, 'Color', 'black', 'HorizontalAlignment', 'center');
text(7, 6, "C", 'FontSize', 12, 'Color', 'black', 'HorizontalAlignment', 'center');
text(5, 10, "E", 'FontSize', 12, 'Color', 'black', 'HorizontalAlignment', 'center');
text(2, 12.5, "F", 'FontSize', 12, 'Color', 'black', 'HorizontalAlignment', 'center');
text(9, 13, "G", 'FontSize', 12, 'Color', 'black', 'HorizontalAlignment', 'center');
text(5, 16, "H", 'FontSize', 12, 'Color', 'black', 'HorizontalAlignment', 'center');

hold off

function slope = slope(stress)

    if stress < 5
        slope = 1;
        return;
    else
        slope = 1 - (stress - 5) * 0.05;
    end

end
