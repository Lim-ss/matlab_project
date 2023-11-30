%{
    f(x)定义在[-1,1]上
    求用节点处二阶导数M表示的三次样条插值函数，令S'(-1) = S'(1) = 0
    采用等距插值
    最终得到等距测试点v上的逐点误差，画出散点图
%} 

%解各内部插值点二次导数值满足的线性方程组A * M = F
n = 2^4;%区间数
h = 2/n;%步长
x = -1:h:1;
y = sin(4 * (x .^ 2)) + (sin(4 * x) .^ 2);
f = @(x)sin(4 * (x .^ 2)) + (sin(4 * x) .^ 2);

A = eye(n+1) * 2;
%M = zeros(n+1,1);%不用预分配，等会会直接创建
F = zeros(n+1,1);

A(1,1) = 2;
A(1,2) = 1;
for i = 2:n
    A(i,i-1) = 0.5;
    A(i,i+1) = 0.5;
end
A(n+1,n) = 1;
A(n+1,n+1) = 2;

F(1) = 6 * (y(2) - y(1)) / h^2; 
for i = 2:n
    F(i) = 3 * (y(i+1) - (2 * y(i)) + y(i-1)) / (h^2);
end
F(n+1) = -6 * (y(n+1) - y(n)) / h^2;

M = A\F;
%{
这段代码画出原函数和插值函数的对比图，用于检查是否插值成功

fplot(f,[-1 1],'r')
hold on
for i = 2:n+1
    S = @(s)(1 / (6 * h)) * ((x(i) - s).^3 * M(i-1) + (s - x(i-1)).^3 * M(i) ...
     + (6 * y(i-1) - M(i-1) * h^2) * (x(i) - s) + (6 * y(i) - M(i) * h^2) * (s - x(i-1)));
    fplot(S,[x(i-1) x(i)],'b')
end
hold off

%} 

%S{1,n+1} = @(z)0;%不知道为什么这样子预分配元胞数组的话，后面会出错
S = cell(1,n+1);%元胞数组，用来存每段多项式函数句柄（注意S(1)没用到，这是因为matlab奇怪的从1开始的下标导致的）
for i = 2:n+1
    S{i} = @(s)(1 / (6 * h)) * ((x(i) - s).^3 * M(i-1) + (s - x(i-1)).^3 * M(i) ...
     + (6 * y(i-1) - M(i-1) * h^2) * (x(i) - s) + (6 * y(i) - M(i) * h^2) * (s - x(i-1)));
end

d = 2000;%等距节点数
v = linspace(-1,1,d);%等距节点
error = zeros(1,d);
error(d) = 0;
segment = 2;%用于指示现在用到第几个S了
for i = 1:d
    if v(i) <= x(segment)
        error(i) = S{segment}(v(i)) - f(v(i));
    else
        segment = segment + 1;
        error(i) = S{segment}(v(i)) - f(v(i));
    end
end

%scatter(v, error);%逐点误差线性散点图
semilogy(v,abs(error));%逐点误差半对数图散点图
