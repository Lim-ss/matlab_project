% 《应力波基础讲义》
% 递减硬化材料中弹塑性波的相互作用（施加的速度大于屈服极限）

v1 = -20;
v2 = 25;
yield_stress = 10;
template = struct('x', 0, 'y', 0);

% 绘制v-sigma图
subplot(1, 4, [1, 2]);
hold on

% 绘制弹性波 aa'  aa"
plot([0, -10], [0, 10], "black", 'LineWidth', 1);
plot([0, 10], [0, 10], "black", 'LineWidth', 1);
% 绘制a"a'''a'a区域
matrix_1 = repmat(template, 11, 16);
matrix_1(1, 1).x = 10;
matrix_1(1, 1).y = 10;

for j = 2:16
    matrix_1(1, j).x = matrix_1(1, j - 1).x + 1;
    matrix_1(1, j).y = matrix_1(1, j - 1).y + slope(matrix_1(1, j - 1).y);
    plot([matrix_1(1, j - 1).x, matrix_1(1, j).x], [matrix_1(1, j - 1).y, matrix_1(1, j).y], "blue", 'LineWidth', 1);
end

for i = 2:11

    matrix_1(i, 1).x = matrix_1(i - 1, 1).x - 1;
    matrix_1(i, 1).y = matrix_1(i - 1, 1).y + slope(matrix_1(i - 1, 1).y);
    plot([matrix_1(i - 1, 1).x, matrix_1(i, 1).x], [matrix_1(i - 1, 1).y, matrix_1(i, 1).y], "blue", 'LineWidth', 1);

    for j = 2:16
        matrix_1(i, j).x = matrix_1(i, j - 1).x + 1;
        matrix_1(i, j).y = matrix_1(i, j - 1).y + slope(matrix_1(i, j - 1).y);
        plot([matrix_1(i, j - 1).x, matrix_1(i, j).x], [matrix_1(i, j - 1).y, matrix_1(i, j).y], "blue", 'LineWidth', 1);
        plot([matrix_1(i - 1, j).x, matrix_1(i, j).x], [matrix_1(i - 1, j).y, matrix_1(i, j).y], "blue", 'LineWidth', 1);
    end

end

% 绘制a'bb'a'''区域

matrix_2 = repmat(template, 11, 11);
matrix_2(1, 1).x = -10;
matrix_2(1, 1).y = 10;

for j = 2:11
    matrix_2(1, j).x = matrix_2(1, j - 1).x + 1;
    matrix_2(1, j).y = matrix_2(1, j - 1).y + slope(matrix_2(1, j - 1).y);
    plot([matrix_2(1, j - 1).x, matrix_2(1, j).x], [matrix_2(1, j - 1).y, matrix_2(1, j).y], "green", 'LineWidth', 1);
end

for i = 2:11

    matrix_2(i, 1).x = matrix_2(i - 1, 1).x - 1;
    matrix_2(i, 1).y = matrix_2(i - 1, 1).y + slope(matrix_2(i - 1, 1).y);
    plot([matrix_2(i - 1, 1).x, matrix_2(i, 1).x], [matrix_2(i - 1, 1).y, matrix_2(i, 1).y], "green", 'LineWidth', 1);

    for j = 2:11
        matrix_2(i, j).x = matrix_2(i, j - 1).x + 1;
        matrix_2(i, j).y = matrix_2(i, j - 1).y + slope(matrix_2(i, j - 1).y);
        plot([matrix_2(i, j - 1).x, matrix_2(i, j).x], [matrix_2(i, j - 1).y, matrix_2(i, j).y], "green", 'LineWidth', 1);
        plot([matrix_2(i - 1, j).x, matrix_2(i, j).x], [matrix_2(i - 1, j).y, matrix_2(i, j).y], "green", 'LineWidth', 1);
    end

end

% 绘制a'''b'dc'区域

matrix_3 = repmat(template, 11, 16);

for i = 1:11
    matrix_3(i, 1).x = matrix_2(i, 11).x;
    matrix_3(i, 1).y = matrix_2(i, 11).y;
end

for j = 1:16
    matrix_3(1, j).x = matrix_1(11, j).x;
    matrix_3(1, j).y = matrix_1(11, j).y;
end

for i = 2:11

    for j = 2:16
        matrix_3(i, j).x = matrix_3(i, j - 1).x + 1;
        matrix_3(i, j).y = matrix_3(i, j - 1).y + slope(matrix_3(i, j - 1).y);
        plot([matrix_3(i, j - 1).x, matrix_3(i, j).x], [matrix_3(i, j - 1).y, matrix_3(i, j).y], "red", 'LineWidth', 1);
        plot([matrix_3(i - 1, j).x, matrix_3(i, j).x], [matrix_3(i - 1, j).y, matrix_3(i, j).y], "red", 'LineWidth', 1);
    end

end

% 绘制坐标轴和标注
text(1, 2, "a", 'FontSize', 12, 'Color', 'black', 'HorizontalAlignment', 'center');
text(matrix_1(1, 1).x + 1, matrix_1(1, 1).y - 1, "a''", 'FontSize', 12, 'Color', 'black', 'HorizontalAlignment', 'center');
text(matrix_1(11, 1).x, matrix_1(11, 1).y - 1, "a'''", 'FontSize', 12, 'Color', 'black', 'HorizontalAlignment', 'center');
text(matrix_1(1, 16).x + 1, matrix_1(1, 16).y, "c", 'FontSize', 12, 'Color', 'black', 'HorizontalAlignment', 'center');
text(matrix_1(11, 16).x, matrix_1(11, 16).y + 1, "c'", 'FontSize', 12, 'Color', 'black', 'HorizontalAlignment', 'center');

text(matrix_2(1, 1).x - 1, matrix_1(1, 1).y - 1, "a'", 'FontSize', 12, 'Color', 'black', 'HorizontalAlignment', 'center');
text(matrix_2(11, 1).x - 1, matrix_1(11, 1).y, "b", 'FontSize', 12, 'Color', 'black', 'HorizontalAlignment', 'center');
text(matrix_2(11, 11).x - 1, matrix_1(11, 11).y + 1, "b'", 'FontSize', 12, 'Color', 'black', 'HorizontalAlignment', 'center');

text(matrix_3(11, 16).x, matrix_1(11, 16).y + 2, "d", 'FontSize', 12, 'Color', 'black', 'HorizontalAlignment', 'center');

plot([0, 0], [-1, 22], "black", 'LineWidth', 1);
plot([-22, 27], [0, 0], "black", 'LineWidth', 1);
plot([-20, -20], [0, matrix_1(11, 1).y], "black--", 'LineWidth', 1);
plot([25, 25], [0, matrix_1(1, 16).y], "black--", 'LineWidth', 1);

text(-20, -1, "v1*", 'FontSize', 12, 'Color', 'black', 'HorizontalAlignment', 'center');
text(25, -1, "v2*", 'FontSize', 12, 'Color', 'black', 'HorizontalAlignment', 'center');

hold off

% 绘制X-t图
subplot(1, 4, [3,4]);
hold on

plot([0, 0], [0, 6], "black", 'LineWidth', 1);
plot([0, 12], [0, 0], "black", 'LineWidth', 1);
plot([10, 10], [0,6], "black--", 'LineWidth', 1);

plot([0, 5], [0, 1.5], "green", 'LineWidth', 1);
plot([10, 5], [0, 1.5], "green", 'LineWidth', 1);

plot([5, 7.5], [1.5, 4], "red", 'LineWidth', 1);
plot([5, 8], [1.5, 4], "red", 'LineWidth', 1);
plot([5, 8.5], [1.5, 4], "red", 'LineWidth', 1);
plot([5, 9], [1.5, 4], "red", 'LineWidth', 1);
plot([5, 9.5], [1.5, 4], "red", 'LineWidth', 1);

plot([5, 2.5], [1.5, 4], "blue", 'LineWidth', 1);
plot([5, 2], [1.5, 4], "blue", 'LineWidth', 1);
plot([5, 1.5], [1.5, 4], "blue", 'LineWidth', 1);
plot([5, 1], [1.5, 4], "blue", 'LineWidth', 1);
plot([5, 0.5], [1.5, 4], "blue", 'LineWidth', 1);

plot([0, 5.5], [0, 5], "yellow", 'LineWidth', 1);
plot([0, 6], [0, 5], "yellow", 'LineWidth', 1);
plot([0, 6.5], [0, 5], "yellow", 'LineWidth', 1);
plot([0, 7], [0, 5], "yellow", 'LineWidth', 1);
plot([0, 7.5], [0, 5], "yellow", 'LineWidth', 1);

plot([10, 4.5], [0, 5], "cyan", 'LineWidth', 1);
plot([10, 4], [0, 5], "cyan", 'LineWidth', 1);
plot([10, 3.5], [0, 5], "cyan", 'LineWidth', 1);
plot([10, 3], [0, 5], "cyan", 'LineWidth', 1);
plot([10, 2.5], [0, 5], "cyan", 'LineWidth', 1);

text(5, 1, "a", 'FontSize', 12, 'Color', 'black', 'HorizontalAlignment', 'center');
text(5, 2.5, "a'''", 'FontSize', 12, 'Color', 'black', 'HorizontalAlignment', 'center');
text(5, 5, "d", 'FontSize', 12, 'Color', 'black', 'HorizontalAlignment', 'center');
text(3, 1.2, "a'", 'FontSize', 12, 'Color', 'black', 'HorizontalAlignment', 'center');
text(7, 1.2, "a''", 'FontSize', 12, 'Color', 'black', 'HorizontalAlignment', 'center');
text(2, 2.5, "b", 'FontSize', 12, 'Color', 'black', 'HorizontalAlignment', 'center');
text(8, 2.5, "c", 'FontSize', 12, 'Color', 'black', 'HorizontalAlignment', 'center');
text(3.7, 3.7, "b'", 'FontSize', 12, 'Color', 'black', 'HorizontalAlignment', 'center');
text(6.3, 3.7, "c'", 'FontSize', 12, 'Color', 'black', 'HorizontalAlignment', 'center');

function slope = slope(stress)

    if stress < 10
        slope = 1;
        return;
    else
        slope = 1 - (stress - 10) * 0.1;
    end

end
