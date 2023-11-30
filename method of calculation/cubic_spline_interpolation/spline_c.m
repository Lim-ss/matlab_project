%{
    在（b）的基础上，令插值区间n不断增大，画出最大误差随n的变化
%} 

%这俩不能放循环内，否则每次循环都会把之前存的值清空
error_max = zeros(1,10);%对应不同n=2^k下的最大误差
x_axis = zeros(1,10);%用于存放画图用的横坐标

for k = 4:10
    n = 2 ^ k;%区间数
    h = 2 / n; %步长
    x = -1:h:1;
    y = sin(4 * (x .^ 2)) + (sin(4 * x) .^ 2);
    f = @(x)sin(4 * (x .^ 2)) + (sin(4 * x) .^ 2);

    A = eye(n + 1) * 2;
    %M = zeros(n+1,1);%不用预分配，等会会直接创建
    F = zeros(n + 1, 1);

    A(1, 1) = 2;
    A(1, 2) = 1;

    for i = 2:n
        A(i, i - 1) = 0.5;
        A(i, i + 1) = 0.5;
    end

    A(n + 1, n) = 1;
    A(n + 1, n + 1) = 2;

    F(1) = 6 * (y(2) - y(1)) / h ^ 2;

    for i = 2:n
        F(i) = 3 * (y(i + 1) - (2 * y(i)) + y(i - 1)) / (h ^ 2);
    end

    F(n + 1) = -6 * (y(n + 1) - y(n)) / h ^ 2;

    M = A \ F;

    S = cell(1, n + 1); %元胞数组，用来存每段多项式函数句柄

    for i = 2:n + 1
        S{i} = @(s)(1 / (6 * h)) * ((x(i) - s) .^ 3 * M(i - 1) + (s - x(i - 1)) .^ 3 * M(i) ...
            + (6 * y(i - 1) - M(i - 1) * h ^ 2) * (x(i) - s) + (6 * y(i) - M(i) * h ^ 2) * (s - x(i - 1)));
    end

    d = 2000; %等距节点数
    v = linspace(-1, 1, d); %等距节点
    error = zeros(1, d);
    error(d) = 0;
    segment = 2; %用于指示现在用到第几个S了
    for i = 1:d

        if v(i) <= x(segment)
            error(i) = S{segment}(v(i)) - f(v(i));
        else
            segment = segment + 1;
            error(i) = S{segment}(v(i)) - f(v(i));
        end
        if error(i) > error_max(k)
            error_max(k) = error(i);
        end
    end

    %scatter(v, error);%逐点误差线性散点图
    %semilogy(v, abs(error)); %逐点误差半对数图散点图

    x_axis(k) = 2^k;
end
%scatter(x_axis(4:10), error_max(4:10));
loglog(x_axis(4:10), error_max(4:10));
