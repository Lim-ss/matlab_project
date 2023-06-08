%{
    使用线性多步法求初值问题y'=xe^(-4x)-4y         y(0) = 0
    采用公式格式为y(n+1) = y(n-1) + [1/3f(n+1) + 4/3f(n) + 1/3f(n-1)]*h       二步方法，只需要再额外一个出发值
    先由二阶的Runge_Kutta法求出出发值
    虽然是隐式公式，但是由于f比较特殊可以改为显式公式y(i+1) = (y(i-1) + h/3*(x(i+1)*exp(-4*x(i+1))) + 4*h/3*f(x(i),y(i)) + h/3*f(x(i-1),y(i-1)))/(1 + 4*h/3)
%} 
h = 0.05;
f = @(a,b) a * exp(-4*a) - 4 * b;
x = 0:h:2;%注意x(1)其实是x0
y(41) = 0;%预分配内存
y(1) = 0;
y(2) = 0.00109396;%由R-K得到的出发值
%y(3) = 0.00935703;
for i = 2:40
    y(i+1) = (y(i-1) + h/3*(x(i+1)*exp(-4*x(i+1))) + 4*h/3*f(x(i),y(i)) + h/3*f(x(i-1),y(i-1)))/(1 + 4*h/3);
end

y_real_value = 0.5 * x.^2 .* exp(-4*x);%精确解，可用来对比

%这段用来对比精确解函数图像和解函数图像
scatter(x, y);
hold on
scatter(x, y_real_value);
hold off

%这段用loglog来判断阶数
error = abs(y - y_real_value);
loglog(x,error);