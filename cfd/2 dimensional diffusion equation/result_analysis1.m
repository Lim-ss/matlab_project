% 用于分析adi方法计算的结果

% 定义数值误差 epsilon(t) = ||u(t) - u'(t)|| 即数值解和解析解的最大误差
% 求epsilon(2)和网格间距的关系（使用对数坐标系）即log(epsilon) ~ log(h)
% 为了方便去log，直接定义网格间距为1/2、1/4、1/8、1/16... 以底数为2取log后，变成-1、-2、-3、-4....

n = 10;
t = 2;
x_plot = zeros(1,n);
y_plot = zeros(1,n);

for i = 1 : n
    delta_x = 1 / 2^i;
    u = analytical_solution(delta_x, t);
    u_ = adi(delta_x, t);
    max = 0;
    for j = 1 : 2^i + 1
        for k = 1 : 2^i + 1
            if max < abs(u(j, k) - u_(j, k))
                max = abs(u(j, k) - u_(j, k));
            end
        end
    end
    x_plot(i) = -i;
    y_plot(i) = log(max) / log(2);
end

plot(x_plot, y_plot);
xlabel('log_2 \Delta_x');
ylabel('log_2 \epsilon(2)');