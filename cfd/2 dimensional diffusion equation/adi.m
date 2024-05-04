% 二维扩散方程 ut = sigma * (uxx + uyy)
% 交替方向隐式方法（Alternating Diretion Implicit, ADI）

function u = adi(delta_x, t_end)
    % 初始化
    Lx = 1;
    Ly = Lx;
    sigma = 1;
    delta_t = 0.01; % 时间步长
    t = 0;
    n = round(Lx / delta_x, 0) + 1; % 网格数四舍五入到小数点后0位，取整
    lambda = sigma * delta_t / (2 * delta_x ^ 2); % 仅用于简化表达式，无实际含义
    %delta_x = input("输入网格间距(x和y的相同)\n"); % 空间步长
    %delta_y = delta_x;
    %t_end = input("输入最终时刻\n");

    % 初始条件（包含初始边界条件）
    u = zeros(n);

    for i = 1:n

        for j = 1:n
            x = ((i - 1) / (n - 1)) * Lx;
            y = ((j - 1) / (n - 1)) * Ly;
            u(i, j) = 20 + 80 * (y - exp(-0.5 * sigma * pi ^ 2 * t) * sin(pi / 2 * x) * sin(pi / 2 * y));
        end

    end

    u_ = zeros(n); % 表示u*
    u_t = zeros(n); % 用于暂时存储t + delta_t时刻u的边界值，在计算u*（即u_）的边界时会用到

    while t < t_end -1e-6
        % 计算下一时刻的边界，结果暂时不存入u而是存入u_t中
        for j = 1:n
            y = ((j - 1) / (n - 1)) * Ly;
            u_t(1, j) = boundaryLeft(y, t + delta_t, sigma);
            u_t(n, j) = boundaryRight(y, t + delta_t, sigma);
        end

        for i = 1:n
            x = ((i - 1) / (n - 1)) * Lx;
            u_t(i, 1) = boundaryDown(x, t + delta_t, sigma);
            u_t(i, n) = boundaryUp(x, t + delta_t, sigma);
        end

        % 计算左右边界（不包括四个角）的u*边界条件（上下边界包括四个角求不出来，也用不到）
        for j = 2:n - 1
            u_(1, j) = 0.5 * ((1 - 2 * lambda) * u(1, j) + (1 + 2 * lambda) * u_t(1, j) + ...
                lambda * (u(1, j - 1) + u(1, j + 1) - u_t(1, j - 1) - u_t(1, j + 1)));
            u_(n, j) = 0.5 * ((1 - 2 * lambda) * u(n, j) + (1 + 2 * lambda) * u_t(n, j) + ...
                lambda * (u(n, j - 1) + u(n, j + 1) - u_t(n, j - 1) - u_t(n, j + 1)));
        end

        % 计算内部区域的u*(u_)
        for j = 2:n - 1
            d = repmat(1 + 2 * lambda, n, 1); % 主对角线上的元素
            e = repmat(-lambda, n, 1); % 上对角线上的元素
            f = repmat(-lambda, n, 1); % 下对角线上的元素
            d(1) = 1;
            d(n) = 1;
            e(1) = 0;
            f(n) = 0;
            b = zeros(n, 1); % 右侧向量
            b(1) = u_(1, j);
            b(n) = u_(n, j);

            for i = 2:n - 1
                b(i) = (1 - 2 * lambda) * u(i, j) + lambda * u(i, j + 1) + lambda * u(i, j - 1);
            end

            u_(:, j) = trisolve(f, d, e, b); % 解三对角方程，返回的结果直接存回u_里
            % 虽然返回的结果是*u1j、u*2j...u*nj 在几何坐标系中是一行，但是在矩阵表示中是一列
        end

        % 计算下一时刻边界的 u （其实已经算过了，直接从u_t里复制过来就行）
        for j = 1:n
            u(1, j) = u_t(1, j);
            u(n, j) = u_t(n, j);
        end

        for i = 1:n
            u(i, 1) = u_t(i, 1);
            u(i, n) = u_t(i, n);
        end

        % 计算下一时刻内部的 u
        for i = 2:n - 1
            d = repmat(1 + 2 * lambda, n, 1); % 主对角线上的元素
            e = repmat(-lambda, n, 1); % 上对角线上的元素
            f = repmat(-lambda, n, 1); % 下对角线上的元素
            d(1) = 1;
            d(n) = 1;
            e(1) = 0;
            f(n) = 0;
            b = zeros(n, 1); % 右侧向量
            b(1) = u(i, 1);
            b(n) = u(i, n);

            for j = 2:n - 1
                b(j) = (1 - 2 * lambda) * u_(i, j) + lambda * u_(i - 1, j) + lambda * u_(i + 1, j);
            end

            u(i, :) = (trisolve(f, d, e, b))'; % 解三对角方程，返回的结果直接存回u里
            % 虽然返回的结果是ui1、ui2...uin 在几何坐标系中是一列，但是在矩阵表示中是一行
        end

        % 更新当前时刻
        t = t + delta_t;
    end

end

% u 的边界条件：
function u = boundaryLeft(y, ~, ~)
    u = 20 + 80 * y;
end

function u = boundaryRight(y, t, sigma)
    u = 20 + 80 * (y - exp(-0.5 * sigma * pi ^ 2 * t) * sin(pi / 2 * y));
end

function u = boundaryDown(~, ~, ~)
    u = 20;
end

function u = boundaryUp(x, t, sigma)
    u = 20 + 80 * (1 - exp(-0.5 * sigma * pi ^ 2 * t) * sin(pi / 2 * x));
end

% 追赶法 LU分解 解三对角方程组，其中abcf分别为下对角线、主对角线、上对角线和右侧向量
% abcf的长度均需要为n，即使a和c本质维度是n-1，a[0]和c[n]不作使用
function x = trisolve(a, b, c, f)
    % 创建辅助数组和结果向量
    n = length(a);
    Beta = zeros(1, n);
    Gamma = zeros(1, n);
    Delta = zeros(1, n);
    y = zeros(1, n);
    x = zeros(n, 1);

    % 三角分解:Ax=f ---> (LU)x = f
    Gamma(1) = b(1);
    Delta(1) = c(1) / Gamma(1);

    for i = 2:n

        Beta(i) = a(i);
        Gamma(i) = b(i) - Beta(i) * Delta(i - 1);

        if i ~= n
            Delta(i) = c(i) / Gamma(i);
        end

    end

    % Ly = f
    y(1) = f(1) / Gamma(1);

    for i = 2:n
        y(i) = (f(i) - Beta(i) * y(i - 1)) / Gamma(i);
    end

    % Ux = y
    x(n) = y(n);

    for i = n - 1:-1:1
        x(i) = y(i) - Delta(i) * x(i + 1);
    end

end
