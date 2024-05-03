% 本题的解析解，用于对比

% 参数
Lx = 1;
Ly = Lx;
sigma = 1; 
delta_x = input("输入网格间距(x和y的相同)\n"); % 空间步长
delta_y = delta_x;
t = input("输入最终时刻\n");
n = round(Lx / delta_x, 0) + 1; % 网格数四舍五入到小数点后0位，取整

% 网格上的解析解
u = zeros(n); 
for i = 1 : n
    for j = 1 : n
        x = ((i - 1) / (n - 1)) * Lx;
        y = ((j - 1) / (n - 1)) * Ly;
        u(i, j) = 20 + 80 * (y - exp(-0.5 * sigma * pi^2 * t) * sin(pi/2 * x) * sin(pi/2 * y));
    end
end