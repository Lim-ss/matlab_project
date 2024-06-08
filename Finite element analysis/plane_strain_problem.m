clear;
clc;

% 读取节点和单元（1-base数组）
[Nodes, Elements] = MeshInput('Nodes.txt', 'Elements.txt');
NodesNum = size(Nodes, 1); % (Vertex)
ElementsNum = size(Elements, 1); % (Triangle)

% 组装总体刚度矩阵
K = zeros(2 * NodesNum, 2 * NodesNum);

for elementIndex = 1:ElementsNum
    i = Elements(elementIndex, 1);
    j = Elements(elementIndex, 2);
    k = Elements(elementIndex, 3);
    Ke = GetElementStiffnessMatrix(elementIndex, Nodes, Elements);

    K(2 * i - 1:2 * i, 2 * i - 1:2 * i) = K(2 * i - 1:2 * i, 2 * i - 1:2 * i) + Ke(1:2, 1:2);
    K(2 * i - 1:2 * i, 2 * j - 1:2 * j) = K(2 * i - 1:2 * i, 2 * j - 1:2 * j) + Ke(1:2, 3:4);
    K(2 * i - 1:2 * i, 2 * k - 1:2 * k) = K(2 * i - 1:2 * i, 2 * k - 1:2 * k) + Ke(1:2, 5:6);

    K(2 * j - 1:2 * j, 2 * i - 1:2 * i) = K(2 * j - 1:2 * j, 2 * i - 1:2 * i) + Ke(3:4, 1:2);
    K(2 * j - 1:2 * j, 2 * j - 1:2 * j) = K(2 * j - 1:2 * j, 2 * j - 1:2 * j) + Ke(3:4, 3:4);
    K(2 * j - 1:2 * j, 2 * k - 1:2 * k) = K(2 * j - 1:2 * j, 2 * k - 1:2 * k) + Ke(3:4, 5:6);

    K(2 * k - 1:2 * k, 2 * i - 1:2 * i) = K(2 * k - 1:2 * k, 2 * i - 1:2 * i) + Ke(5:6, 1:2);
    K(2 * k - 1:2 * k, 2 * j - 1:2 * j) = K(2 * k - 1:2 * k, 2 * j - 1:2 * j) + Ke(5:6, 3:4);
    K(2 * k - 1:2 * k, 2 * k - 1:2 * k) = K(2 * k - 1:2 * k, 2 * k - 1:2 * k) + Ke(5:6, 5:6);

end

clear i j k Ke elementIndex;

% 获取边界
upline = GetLine(Nodes, Elements, 0, 3.3, 5.196, 5.196);
downline = GetLine(Nodes, Elements, 3, 6, 0, 0);
rightline = GetLine(Nodes, Elements, 3.2, 6, 5.196, 0);
leftline = GetLine(Nodes, Elements, 0, 0, 3, 5.196);
arc = GetCircleArc(Nodes, 0, 0, 3);

% 施加力和位移边界条件
F = zeros(2 * NodesNum, 1);
F = FixedStressBoundaryConditions(upline, Nodes, F, 0, 0, 0, 5000);
[K, F] = FixedDisplacementBoundaryConditions(leftline, Nodes, K, F, 0, "x");
[K, F] = FixedDisplacementBoundaryConditions(downline, Nodes, K, F, 0, "x");
[K, F] = FixedDisplacementBoundaryConditions(downline, Nodes, K, F, 0, "y");

% 求解方程组 Kd=F
d = K \ F; % 注意不要写成 F \ K

% 绘制网格
hold on;
PlotMesh(Nodes, Elements, "black");
plotLine(upline, Nodes, "red");
plotLine(downline, Nodes, "blue");
plotLine(rightline, Nodes, "green");
plotLine(leftline, Nodes, "blue");
plotCircleArc(arc, Nodes, "red");
%plotIndex(Nodes);
axis equal;

newNodes = GetNewPosition(Nodes, d);
PlotMesh(newNodes, Elements, "red");
hold off;

% 获得单元刚度矩阵Ke
function Ke = GetElementStiffnessMatrix(elementIndex, Nodes, Elements)
    % 选项
    option1 = 2; % 平面应力问题 —— 1，平面应变问题 —— 2
    % 常数
    E = 2E11;
    nu = 0.3;

    xi = Nodes(Elements(elementIndex, 1), 1);
    yi = Nodes(Elements(elementIndex, 1), 2);
    xj = Nodes(Elements(elementIndex, 2), 1);
    yj = Nodes(Elements(elementIndex, 2), 2);
    xk = Nodes(Elements(elementIndex, 3), 1);
    yk = Nodes(Elements(elementIndex, 3), 2);

    cross_product = (xj - xi) * (yk - yi) - (yj - yi) * (xk - xi);

    if cross_product > 0
        type = "counter-clockwise";
    else
        type = "clockwise";
    end

    if (type == "clockwise")
        t = xi;
        xi = xj;
        xj = t;
        t = yi;
        yi = yj;
        yj = t;
    end

    t = 1;
    bi = yj - yk;
    bj = yk - yi;
    bk = yi - yj;
    ci = xk - xj;
    cj = xi - xk;
    ck = xj - xi;

    delta = abs(((xj - xi) * (yk - yi) - (xk - xi) * (yj - yi)) / 2);

    B = [bi, 0, bj, 0, bk, 0;
         0, ci, 0, cj, 0, ck;
         ci, bi, cj, bj, ck, bk] / (2 * delta);

    if option1 == 1
        % 平面应力问题
        D = E / (1 - nu ^ 2) * [1, nu, 0; nu, 1, 0; 0, 0, (1 - nu) / 2];
    elseif option1 == 2
        E_ = E / (1 - nu ^ 2);
        nu_ = nu / (1 - nu);
        D = E_ / (1 - nu_ ^ 2) * [1, nu_, 0; nu_, 1, 0; 0, 0, (1 - nu_) / 2];
    else
        disp("刚度矩阵非法选项")
    end

    Ke = t * delta * B' * D * B;

    if (type == "clockwise")
        % 交换13列和24列
        T = Ke(1:6, 1:2);
        Ke(1:6, 1:2) = Ke(1:6, 3:4);
        Ke(1:6, 3:4) = T;
    end

    return;
end

function PlotMesh(Nodes, Elements, color)
    ElementsNum = size(Elements, 1);

    for elementIndex = 1:ElementsNum
        xi = Nodes(Elements(elementIndex, 1), 1);
        yi = Nodes(Elements(elementIndex, 1), 2);
        xj = Nodes(Elements(elementIndex, 2), 1);
        yj = Nodes(Elements(elementIndex, 2), 2);
        xk = Nodes(Elements(elementIndex, 3), 1);
        yk = Nodes(Elements(elementIndex, 3), 2);
        plot([xi, xj], [yi, yj], color, 'LineWidth', 1);
        plot([xj, xk], [yj, yk], color, 'LineWidth', 1);
        plot([xk, xi], [yk, yi], color, 'LineWidth', 1);
    end

end

function Line = GetLine(Nodes, Elements, x1, x2, y1, y2)

    LineNodes = zeros(1, 1); % 避免警告
    NodesNum = size(Nodes, 1);
    num = 0; % 线段上含有的点的个数
    tolerance = 1e-5;

    for nodeIndex = 1:NodesNum
        x = Nodes(nodeIndex, 1);
        y = Nodes(nodeIndex, 2);

        if (x == x1 && y == y1) || (x == x2 && y == y2)
            num = num + 1;
            LineNodes(num) = nodeIndex;
            continue;
        end

        if x1 == x2

            if x == x1 && ((y1 < y && y < y2) || (y2 < y && y < y1))
                num = num + 1;
                LineNodes(num) = nodeIndex;
                continue;
            end

        end

        if y1 == y2

            if y == y1 && ((x1 < x && x < x2) || (x2 < x && x < x1))
                num = num + 1;
                LineNodes(num) = nodeIndex;
                continue;
            end

        end

        if abs((x -x1) / (x2 - x1) - (y - y1) / (y2 - y1)) < tolerance
            num = num + 1;
            LineNodes(num) = nodeIndex;
            continue;
        end

    end

    % 以x为第一关键字y为第二关键字进行排序
    if x1 ~= x2

        for i = 1:num - 1

            for j = 1:num - i

                if Nodes(LineNodes(j), 1) > Nodes(LineNodes(j + 1), 1)

                    t = LineNodes(j);
                    LineNodes(j) = LineNodes(j + 1);
                    LineNodes(j + 1) = t;
                end

            end

        end

    else

        for i = 1:num - 1

            for j = 1:num - i

                if Nodes(LineNodes(j), 2) > Nodes(LineNodes(j + 1), 2)

                    t = LineNodes(j);
                    LineNodes(j) = LineNodes(j + 1);
                    LineNodes(j + 1) = t;
                end

            end

        end

    end

    %查找并按顺序存放线段上的三角形

    isInLine = zeros(NodesNum);
    rank = ones(NodesNum) * 9999;

    for i = 1:num
        isInLine(LineNodes(i)) = 1;
        rank(LineNodes(i)) = i;
    end

    ElementsNum = size(Elements, 1);
    LineElements = zeros(1, num - 1);
    num2 = 0;

    for i = 1:ElementsNum
        inLineNum = isInLine(Elements(i, 1)) + isInLine(Elements(i, 2)) + isInLine(Elements(i, 3));

        if inLineNum == 2
            minRank = min(rank(Elements(i, 1)), min(rank(Elements(i, 2)), rank(Elements(i, 3))));
            num2 = num2 + 1;
            LineElements(minRank) = i;
        end

    end

    Line = {LineNodes, LineElements};

    return;
end

function Arc = GetCircleArc(Nodes, x1, y1, r)
    Arc = zeros(1, 1); % 避免警告
    angle = zeros(1, 1); % 每个弧上的节点对应的角度，用于排序
    NodesNum = size(Nodes, 1);
    num = 0; % 线段上含有的点的个数
    tolerance = 1e-5;

    for nodeIndex = 1:NodesNum
        x = Nodes(nodeIndex, 1);
        y = Nodes(nodeIndex, 2);

        if (x - x1) ^ 2 + (y - y1) ^ 2 - r ^ 2 < tolerance
            num = num + 1;
            Arc(num) = nodeIndex;
            angle(num) = rad2deg(atan2(y - y1, x - x1));
        end

    end

    % 按照角度进行排序
    for i = 1:num - 1

        for j = 1:num - i

            if angle(j) > angle(j + 1)
                t = angle(j);
                angle(j) = angle(j + 1);
                angle(j + 1) = t;
                t = Arc(j);
                Arc(j) = Arc(j + 1);
                Arc(j + 1) = t;
            end

        end

    end

    return;
end

function plotLine(Line, Nodes, color)

    LineNodes = Line{1};
    num = size(LineNodes, 2);
    x1 = Nodes(LineNodes(1), 1);
    x2 = Nodes(LineNodes(num), 1);
    y1 = Nodes(LineNodes(1), 2);
    y2 = Nodes(LineNodes(num), 2);
    plot([x1, x2], [y1, y2], color, 'LineWidth', 1);

end

function plotCircleArc(Arc, Nodes, color)
    num = size(Arc, 2);

    for i = 1:num - 1
        x1 = Nodes(Arc(i), 1);
        y1 = Nodes(Arc(i), 2);
        x2 = Nodes(Arc(i + 1), 1);
        y2 = Nodes(Arc(i + 1), 2);

        plot([x1, x2], [y1, y2], color, 'LineWidth', 1);
    end

end

function [Nodes, Elements] = MeshInput(fnodes, felements)
    [f, message] = fopen(fnodes, 'r');

    if f == -1
        disp(message);
    end

    C2 = textscan(f, '%f %f');
    Nodes(:, 1) = C2 {1};
    Nodes(:, 2) = C2 {2};
    fclose(f);

    [f, message] = fopen(felements, 'r');

    if f == -1
        disp(message);
    end

    C3 = textscan(f, '%f %f %f');
    Elements(:, 1) = C3 {1};
    Elements(:, 2) = C3 {2};
    Elements(:, 3) = C3 {3};
    fclose(f);

    return;
end

function F = FixedStressBoundaryConditions(Line, Nodes, F, stress1X, stress2X, stress1Y, stress2Y)

    % stress1和stress2优先根据x进行先后对应
    LineNodes = Line{1};
    LineElements = Line{2};
    num = size(LineNodes, 2);

    x1 = Nodes(LineNodes(1), 1);
    x2 = Nodes(LineNodes(num), 1);
    y1 = Nodes(LineNodes(1), 2);
    y2 = Nodes(LineNodes(num), 2);

    if (x1 ~= x2)

        slopeX = (stress2X - stress1X) / (x2 - x1);
        slopeY = (stress2Y - stress1Y) / (x2 - x1);

        for i = 1:num - 1
            xa = Nodes(LineNodes(i), 1);
            xb = Nodes(LineNodes(i + 1), 1);
            stressX = (stress1X + slopeX * (xa - x1) + stress1X + slopeX * (xb - x1)) * (xb - xa) / 2;
            stressY = (stress1Y + slopeY * (xa - x1) + stress1Y + slopeY * (xb - x1)) * (xb - xa) / 2;

            % n1 = LineElements(i);
            % n2 = LineElements(i);
            % n3 = LineElements(i);

            % F(n1 * 2 - 1) = F(n1 * 2 - 1) + stressX / 3;
            % F(n2 * 2 - 1) = F(n2 * 2 - 1) + stressX / 3;
            % F(n3 * 2 - 1) = F(n3 * 2 - 1) + stressX / 3;

            % F(n1 * 2) = F(n1 * 2) + stressY / 3;
            % F(n2 * 2) = F(n2 * 2) + stressY / 3;
            % F(n3 * 2) = F(n3 * 2) + stressY / 3;

            F(LineNodes(i) * 2 - 1) = F(LineNodes(i) * 2 - 1) + stressX / 2;
            F(LineNodes(i + 1) * 2 - 1) = F(LineNodes(i + 1) * 2 - 1) + stressX / 2;
            F(LineNodes(i) * 2) = F(LineNodes(i) * 2) + stressY / 2;
            F(LineNodes(i + 1) * 2) = F(LineNodes(i + 1) * 2) + stressY / 2;
        end

    else
        slopeX = (stress2X - stress1X) / (y2 - y1);
        slopeY = (stress2Y - stress1Y) / (y2 - y1);

        for i = 1:num - 1
            ya = Nodes(LineNodes(i), 2);
            yb = Nodes(LineNodes(i + 1), 2);
            stressX = (stress1X + slopeX * (ya - y1) + stress1X + slopeX * (yb - y1)) * (yb - ya) / 2;
            stressY = (stress1Y + slopeY * (ya - y1) + stress1Y + slopeY * (yb - y1)) * (yb - ya) / 2;

            % n1 = LineElements(i);
            % n2 = LineElements(i);
            % n3 = LineElements(i);

            % F(n1 * 2 - 1) = F(n1 * 2 - 1) + stressX / 3;
            % F(n2 * 2 - 1) = F(n2 * 2 - 1) + stressX / 3;
            % F(n3 * 2 - 1) = F(n3 * 2 - 1) + stressX / 3;

            % F(n1 * 2) = F(n1 * 2) + stressY / 3;
            % F(n2 * 2) = F(n2 * 2) + stressY / 3;
            % F(n3 * 2) = F(n3 * 2) + stressY / 3;

            F(LineNodes(i) * 2 - 1) = F(LineNodes(i) * 2 - 1) + stressX / 2;
            F(LineNodes(i + 1) * 2 - 1) = F(LineNodes(i + 1) * 2 - 1) + stressX / 2;
            F(LineNodes(i) * 2) = F(LineNodes(i) * 2) + stressY / 2;
            F(LineNodes(i + 1) * 2) = F(LineNodes(i + 1) * 2) + stressY / 2;
        end

    end

    return;
end

function [K, F] = FixedDisplacementBoundaryConditions(Line, Nodes, K, F, displacement, option)

    % option = "x"代表固定x方向位移，option = "y"代表固定y方向位移
    LineNodes = Line{1};
    num = size(LineNodes, 2);

    if (option == "x")

        for i = 1:num

            for j = 1:size(Nodes, 1) * 2
                F(j) = F(j) - K(2 * LineNodes(i) - 1, j) * displacement;
                K(j, 2 * LineNodes(i) - 1) = 0;
                K(2 * LineNodes(i) - 1, j) = 0;
            end

            K(2 * LineNodes(i) - 1, 2 * LineNodes(i) - 1) = 1;
            F(2 * LineNodes(i) - 1) = displacement;
        end

    elseif (option == "y")

        for i = 1:num

            for j = 1:size(Nodes, 1) * 2
                F(j) = F(j) - K(2 * LineNodes(i), j) * displacement;
                K(j, 2 * LineNodes(i)) = 0;
                K(2 * LineNodes(i), j) = 0;
            end

            K(2 * LineNodes(i), 2 * LineNodes(i)) = 1;
            F(2 * LineNodes(i)) = displacement;
        end

    else
        disp("固定位移边界条件非法输入");
    end

end

function plotIndex(Nodes)
    NodesNum = size(Nodes, 1);

    for i = 1:NodesNum
        text(Nodes(i, 1), Nodes(i, 2), num2str(i), 'FontSize', 10, 'Color', 'black', 'HorizontalAlignment', 'center');
    end

end

function NewNodes = GetNewPosition(Nodes, d)
    times = 3e6;
    num = size(Nodes, 1);
    NewNodes = zeros(num, 2);

    for i = 1:num
        NewNodes(i, 1) = Nodes(i, 1) + d(2 * i - 1) * times;
        NewNodes(i, 2) = Nodes(i, 2) + d(2 * i) * times;
    end

    return;
end
