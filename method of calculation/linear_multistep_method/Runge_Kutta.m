%{
    二阶R-K方法，用于计算线性多步法的起步值
    y'=xe^(-4x)-4y         y(0) = 0
    定义f = xe^(-4x)-4y
%} 
h = 0.05;%步长
f = @(a,b) a * exp(-4*a) - 4 * b;
x = 0:h:2;%注意x(1)其实是x0
y(41) = 0;
y(1) = 0;
for i = 1:40
    k1 = h * f(x(i),y(i));
    k2 = h * f(x(i)+ (h * 2/3), y(i) + (k1 * 2/3));
    y(i+1) = y(i) + (k1 + 3*k2)/4;
end

%最后得到的y0，y1，y2……按顺序存在y数组中，直接查看即可

y_real_value = 0.5 * x.^2 .* exp(-4*x);%精确解，可用来对比

scatter(x, y);
hold on
scatter(x, y_real_value);
hold off