% 《应力波基础》（第二版）P49页
% 阶梯弹性杆中波的反射和透射

subplot(2, 2, [1, 2]);
hold on;
% 绘制坐标轴
plot([0, 0], [-1, 6], "black");
plot([1, 1], [0, 6], "black--");
plot([2, 2], [0, 6], "black--");
plot([3, 3], [0, 6], "black--");
plot([4, 4], [0, 6], "black--");
plot([5, 5], [0, 6], "black--");
plot([-2, 6], [0, 0], "black");

% 绘制行波
plotWaveLine(0, 0, 1);
plotWaveLine(0, 0, -1);

stress = zeros(6, 6); % (0s、1s、...、5s时刻)各区域内的应力大小
strength = zeros(6, 5); % 应力波的改变强度
transmission_right = zeros(1, 5);
reflection_right = zeros(1, 5);
transmission_left = zeros(1, 5);
reflection_left = zeros(1, 5);

% 计算各处的投射/反射系数
for i = 1:5
    k1 = ((7 - i) / (6 - i)) ^ 2;
    k2 = ((6 - i) / (7 - i)) ^ 2;
    transmission_right(i) = (2 * k1) / (1 + k1);
    reflection_right(i) = (1 - k1) / (1 + k1);
    transmission_left(i) = (2 * k2) / (1 + k2);
    reflection_left(i) = (1 - k2) / (1 + k2);
end

% 计算各行波的强度（改变倍数）
strength(1, 1) = -1 * reflection_right(1);
strength(1, 2) = -1 * transmission_right(1);

for i = 2:5

    for j = 1:6

        if mod(i + j, 2) == 1
            % 右行波
            if j == 1
                strength(i, j) = 0;
            else
                strength(i, j) = strength(i - 1, j - 1) * transmission_right(j - 1) + strength(i - 1, j) * reflection_left(j - 1);
            end

        else
            %左行波
            if j == 6
                strength(i, j) = 0;
            else
                strength(i, j) = strength(i - 1, j + 1) * transmission_left(j) + strength(i - 1, j) * reflection_right(j);
            end

        end

    end

end

% 在图中用数字标记行波的强度
for i = 1:5

    for j = 1:6

        if strength(i, j) ~= 0
            text(j - 1.5, i - 0.5, num2str(strength(i, j)), 'FontSize', 9, 'Color', 'r', 'HorizontalAlignment', 'center');
        end

    end

end

% 计算各处的应力(但是题目要求的不是画这个)
for j = 1:6
    stress(1, j) = 0;
end

for i = 2:6

    for j = 1:6
        stress(i, j) = stress(i - 1, j) + strength(i - 1, j);
    end

end

hold off;

subplot(2, 2, 3);
hold on
% 绘制t=3->t=4的波的强度图
plot([0, 0], [-3, 1], "black");
plot([1, 1], [-3, 1], "black--");
plot([2, 2], [-3, 1], "black--");
plot([3, 3], [-3, 1], "black--");
plot([4, 4], [-3, 1], "black--");
plot([5, 5], [-3, 1], "black--");
plot([0, 6], [0, 0], "black");

for i = 2:6
    plot([i - 2, i - 1], [strength(4, i), strength(4, i)], "black", 'LineWidth', 2);
    text(i - 1.5, strength(4, i) + 0.5, num2str(strength(4, i)), 'FontSize', 9, 'Color', 'r', 'HorizontalAlignment', 'center');
end

plot([0, 0], [0, strength(4, 2)], "black", 'LineWidth', 2);
plot([5, 5], [0, strength(4, 6)], "black", 'LineWidth', 2);

for i = 2:5
    plot([i - 1, i - 1], [strength(4, i), strength(4, i + 1)], "black", 'LineWidth', 2);
end

hold off

subplot(2, 2, 4);
hold on
% 绘制t=4->t=5的波的强度图
plot([0, 0], [-4.5, 2], "black");
plot([1, 1], [-4.5, 2], "black--");
plot([2, 2], [-4.5, 2], "black--");
plot([3, 3], [-4.5, 2], "black--");
plot([4, 4], [-4.5, 2], "black--");
plot([5, 5], [-4.5, 2], "black--");
plot([0, 6], [0, 0], "black");

for i = 2:6
    plot([i - 2, i - 1], [strength(5, i), strength(5, i)], "black", 'LineWidth', 2);
    text(i - 1.5, strength(5, i) + 0.5, num2str(strength(5, i)), 'FontSize', 9, 'Color', 'r', 'HorizontalAlignment', 'center');
end

plot([0, 0], [0, strength(5, 2)], "black", 'LineWidth', 2);
plot([5, 5], [0, strength(5, 6)], "black", 'LineWidth', 2);

for i = 2:5
    plot([i - 1, i - 1], [strength(5, i), strength(5, i + 1)], "black", 'LineWidth', 2);
end

hold off

% 传入一个位置和时间，并指明左行波还是右行波，diretion = 1代表右行，diretion = -1代表左行
function plotWaveLine(x, t, diretion)

    % 设置图标的参数
    left_x = 0;
    right_x = 5;
    time_end = 5;
    k = 1; % 斜率，绘制的时候y轴乘上k即可
    lineColor = "blue";

    if diretion == 1

        if x == right_x
            return;
        end

        if t == time_end
            plot([x, x + 0.5], [k * t, k * (t + 0.5)], lineColor);
        else
            plot([x, x + 1], [k * t, k * (t + 1)], lineColor);
            plotWaveLine(x + 1, t + 1, 1);
            plotWaveLine(x + 1, t + 1, -1);
        end

        return;
    else
        % diretion == -1
        if x == left_x
            plot([x, x - 1], [k * t, k * (t + 1)], lineColor);
            return;
        end

        if t == time_end
            plot([x, x - 0.5], [k * t, k * (t + 0.5)], lineColor);
        else
            plot([x, x - 1], [k * t, k * (t + 1)], lineColor);
            plotWaveLine(x - 1, t + 1, 1);
            plotWaveLine(x - 1, t + 1, -1);
        end

        return;
    end

end
