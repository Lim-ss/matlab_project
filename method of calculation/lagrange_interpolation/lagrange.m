%{
    [0,1]上周期函数基于等距节点的拉格朗日插值
%} 
n = 2^6;%区间数
x = 0:1/n:(n-1)/n;

f = @(a) sin(2*pi*a) .* exp(cos(2*pi*a));
y = f(x);

S = cell(1,n);%元胞数组，用来存放基函数,注意S(i)里面放的是l_(i-1)
for i = 1:n
    k = i - 1;
    if mod(n,2) == 1
        S{i} = @(s)(-1)^k * sin(n*pi*s) .* csc(pi*(s - k/n)) / n;
    else
        S{i} = @(s)(-1)^k * sin(n*pi*s) .* cot(pi*(s - k/n)) / n;
    end
end

d = 1000;%等距节点数
v = linspace(0,1,d);%等距节点
result = zeros(1,d);%存放插值函数算出的结果
error = zeros(1,d);
real_value = zeros(1,d);%保存测试点的真实值，用于与result比较，检查插值是否正确

for i = 1:d
    for j = 1:n
        result(i) = result(i) + S{j}(v(i)) * y(j);
    end
    error(i) = result(i) - f(v(i));
    real_value(i) = f(v(i));
end

%{
    分别画出真值和拟合值的散点图，用于debug时检查正确性
scatter(v, result);
hold on
scatter(v, real_value);
hold off
%}
semilogy(v, abs(error));
