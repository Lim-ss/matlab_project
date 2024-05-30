%《振动理论及应用》（第五版）P123 题4.34
% 利用数值方法绘制无阻尼受迫振动下的响应的图像
% 控制方程为 4x" + 2000x = F(t)
% 迭代公式：x(i+1) = 2x(i) - x(i-1) + h^2 * (1/4 F(t) - 500 x(i))

h = 0.01; % 时间间隔
t_end = 0.5; % 最大计算时间

n = round(t_end / h, 0) + 1;
x = zeros(1, n);
x(1) = 0;
x(2) = x(1) + h * 0 + 0.5 * h ^ 2 * 0.25* F(0);

i = 3;

while (i - 1) * h <= t_end
    t = (i - 1) * h;
    x(i) = 2 * x(i - 1) - x(i - 2) + h ^ 2 * (1/4 * F(t) - 500 * x(i - 1));
    i = i + 1;
end

t = 0 : h : t_end;
plot(t, x, 'r');
grid on;
title('x-t');
xlabel('t');
ylabel('x');


function F = F(t)

    if t <= 0.1
        F = 100;
    elseif t <= 0.2
        F = 200 - 1000 * t;
    else
        F = 0;
    end

end
