% 用于分析adi方法计算的结果

% 定义计算误差 E = ||u^(n+1) - u^(n)|| 即数值解前后时刻的最大差异
% 求log(E)随时间的变换曲线，这里以10为底求对数

delta_x = 0.1;
delta_t = 0.1; % 注意这个delta_t不是adi的时间步长，只用于分析采样
n = round(1 / delta_x, 0) + 1;
m = 2 / delta_t + 1;
x_plot = zeros(1,m);
y_plot = zeros(1,m);

u1 = adi(delta_x, 0);
for i = 1 : m
    t = 0 + delta_t * (i - 1);
    u2 = adi(delta_x, t + delta_t);
    max = 0;
    for j = 1 : n
        for k = 1 : n
            if max < abs(u2(j, k) - u1(j, k))
                max = abs(u2(j, k) - u1(j, k));
            end
        end
    end
    u1 = u2;
    x_plot(i) = t;
    y_plot(i) = log(max)/log(10);
end

plot(x_plot, y_plot);
xlabel('t');
ylabel('logE');