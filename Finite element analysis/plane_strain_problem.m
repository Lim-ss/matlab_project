clear;
clc;

% 读取节点和单元（1-base数组）
[Nodes, Elements] = MeshInput('Nodes2.txt', 'Elements2.txt');
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
upline = GetLine(Nodes, 0, 3.3, 5.196, 5.196);
downline = GetLine(Nodes, 3, 6, 0, 0);
rightline = GetLine(Nodes, 3.2, 6, 5.196, 0);
leftline = GetLine(Nodes, 0, 0, 3, 5.196);
arc = GetCircleArc(Nodes, 0, 0, 3);

% 施加力和位移边界条件
F = zeros(2 * NodesNum, 1);
F = FixedStressBoundaryConditions(upline, Nodes, F, 0, 0, 0, 5000);
[K, F] = FixedDisplacementBoundaryConditions(leftline, Nodes, K, F, 0, "x");
[K, F] = FixedDisplacementBoundaryConditions(downline, Nodes, K, F, 0, "x");
[K, F] = FixedDisplacementBoundaryConditions(downline, Nodes, K, F, 0, "y");

% 求解方程组 Kd=F
% d中元素按照u,v的顺序排列
d = K \ F; % 注意不要写成 F \ K

%重新组织d的格式，让uv在同一行，看起来更清晰
for i = 1:NodesNum
    d(i, 1) = d(2 * i - 1);
    d(i, 2) = d(2 * i);
end

d = d(1:NodesNum, 1:2);
clear i;

[strain, stress] = GetStrainAndStress(Nodes, Elements, d);
[maxStress, belongElement] = GetMaxStress(stress, "Mises");
newNodes = GetNewPosition(Nodes, d);

% 绘制网格/输出信息
hold on;
axis equal;
PlotMesh(Nodes, Elements, "black");
plotLine(upline, Nodes, "red");
plotLine(downline, Nodes, "blue");
plotLine(rightline, Nodes, "green");
plotLine(leftline, Nodes, "blue");
plotCircleArc(arc, Nodes, "red");
%plotNodesIndex(Nodes,"All");
%plotElementsIndex(Nodes, Elements,"All");
PlotMesh(newNodes, Elements, "red");
OutPutRegionData(Nodes, Elements, 3, 1, 0.2);
hold off;

% 获得单元刚度矩阵Ke
function Ke = GetElementStiffnessMatrix(elementIndex, Nodes, Elements)
    % 选项，在GetStrainAndStress()中也需要修改
    option1 = 2; % 平面应力问题 —— 1，平面应变问题 —— 2
    % 常数
    E = 2E11;
    nu = 0.3;
    t = 1;

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

function Line = GetLine(Nodes, x1, x2, y1, y2)

    Line = zeros(1, 1); % 避免警告
    NodesNum = size(Nodes, 1);
    num = 0; % 线段上含有的点的个数
    tolerance = 1e-5;

    for nodeIndex = 1:NodesNum
        x = Nodes(nodeIndex, 1);
        y = Nodes(nodeIndex, 2);

        if (x == x1 && y == y1) || (x == x2 && y == y2)
            num = num + 1;
            Line(num) = nodeIndex;
            continue;
        end

        if x1 == x2

            if x == x1 && ((y1 < y && y < y2) || (y2 < y && y < y1))
                num = num + 1;
                Line(num) = nodeIndex;
                continue;
            end

        end

        if y1 == y2

            if y == y1 && ((x1 < x && x < x2) || (x2 < x && x < x1))
                num = num + 1;
                Line(num) = nodeIndex;
                continue;
            end

        end

        if abs((x -x1) / (x2 - x1) - (y - y1) / (y2 - y1)) < tolerance
            num = num + 1;
            Line(num) = nodeIndex;
            continue;
        end

    end

    % 以x为第一关键字y为第二关键字进行排序
    if x1 ~= x2

        for i = 1:num - 1

            for j = 1:num - i

                if Nodes(Line(j), 1) > Nodes(Line(j + 1), 1)

                    t = Line(j);
                    Line(j) = Line(j + 1);
                    Line(j + 1) = t;
                end

            end

        end

    else

        for i = 1:num - 1

            for j = 1:num - i

                if Nodes(Line(j), 2) > Nodes(Line(j + 1), 2)

                    t = Line(j);
                    Line(j) = Line(j + 1);
                    Line(j + 1) = t;
                end

            end

        end

    end

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

    num = size(Line, 2);
    x1 = Nodes(Line(1), 1);
    x2 = Nodes(Line(num), 1);
    y1 = Nodes(Line(1), 2);
    y2 = Nodes(Line(num), 2);
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
    num = size(Line, 2);

    x1 = Nodes(Line(1), 1);
    x2 = Nodes(Line(num), 1);
    y1 = Nodes(Line(1), 2);
    y2 = Nodes(Line(num), 2);

    if (x1 ~= x2)

        slopeX = (stress2X - stress1X) / (x2 - x1);
        slopeY = (stress2Y - stress1Y) / (x2 - x1);

        for i = 1:num - 1
            xa = Nodes(Line(i), 1);
            xb = Nodes(Line(i + 1), 1);
            ForceX = (stress1X + slopeX * (xa - x1) + stress1X + slopeX * (xb - x1)) * (xb - xa) / 2;
            ForceY = (stress1Y + slopeY * (xa - x1) + stress1Y + slopeY * (xb - x1)) * (xb - xa) / 2;

            F(Line(i) * 2 - 1) = F(Line(i) * 2 - 1) + ForceX / 2;
            F(Line(i + 1) * 2 - 1) = F(Line(i + 1) * 2 - 1) + ForceX / 2;
            F(Line(i) * 2) = F(Line(i) * 2) + ForceY / 2;
            F(Line(i + 1) * 2) = F(Line(i + 1) * 2) + ForceY / 2;
        end

    else
        slopeX = (stress2X - stress1X) / (y2 - y1);
        slopeY = (stress2Y - stress1Y) / (y2 - y1);

        for i = 1:num - 1
            ya = Nodes(Line(i), 2);
            yb = Nodes(Line(i + 1), 2);
            ForceX = (stress1X + slopeX * (ya - y1) + stress1X + slopeX * (yb - y1)) * (yb - ya) / 2;
            ForceY = (stress1Y + slopeY * (ya - y1) + stress1Y + slopeY * (yb - y1)) * (yb - ya) / 2;

            F(Line(i) * 2 - 1) = F(Line(i) * 2 - 1) + ForceX / 2;
            F(Line(i + 1) * 2 - 1) = F(Line(i + 1) * 2 - 1) + ForceX / 2;
            F(Line(i) * 2) = F(Line(i) * 2) + ForceY / 2;
            F(Line(i + 1) * 2) = F(Line(i + 1) * 2) + ForceY / 2;
        end

    end

    return;
end

function [K, F] = FixedDisplacementBoundaryConditions(Line, Nodes, K, F, displacement, option)

    % option = "x"代表固定x方向位移，option = "y"代表固定y方向位移
    num = size(Line, 2);

    if (option == "x")

        for i = 1:num

            for j = 1:size(Nodes, 1) * 2
                F(j) = F(j) - K(2 * Line(i) - 1, j) * displacement;
                K(j, 2 * Line(i) - 1) = 0;
                K(2 * Line(i) - 1, j) = 0;
            end

            K(2 * Line(i) - 1, 2 * Line(i) - 1) = 1;
            F(2 * Line(i) - 1) = displacement;
        end

    elseif (option == "y")

        for i = 1:num

            for j = 1:size(Nodes, 1) * 2
                F(j) = F(j) - K(2 * Line(i), j) * displacement;
                K(j, 2 * Line(i)) = 0;
                K(2 * Line(i), j) = 0;
            end

            K(2 * Line(i), 2 * Line(i)) = 1;
            F(2 * Line(i)) = displacement;
        end

    else
        disp("固定位移边界条件非法输入");
    end

end

function plotNodesIndex(Nodes, option)

    NodesNum = size(Nodes, 1);

    if option == "All"

        for i = 1:NodesNum
            text(Nodes(i, 1), Nodes(i, 2), num2str(i), 'FontSize', 10, 'Color', 'black', 'HorizontalAlignment', 'center');
        end

    elseif option == "region"

        % 只显示特定区域附近的编号
        xx = 2;
        yy = 0;
        radius = 0.2;

        for i = 1:NodesNum
            x = Nodes(i, 1);
            y = Nodes(i, 2);

            if abs(x - xx) < radius && abs(y - yy) < radius
                text(x, y, num2str(i), 'FontSize', 10, 'Color', 'black', 'HorizontalAlignment', 'center');
            end

        end

    else

        i = str2double(option);
        text(Nodes(i, 1), Nodes(i, 2), num2str(i), 'FontSize', 10, 'Color', 'black', 'HorizontalAlignment', 'center');

    end

end

function plotElementsIndex(Nodes, Elements, option)

    ElementsNum = size(Elements, 1);

    if option == "All"

        for i = 1:ElementsNum
            x = (Nodes(Elements(i, 1), 1) + Nodes(Elements(i, 2), 1) + Nodes(Elements(i, 3), 1)) / 3;
            y = (Nodes(Elements(i, 1), 2) + Nodes(Elements(i, 2), 2) + Nodes(Elements(i, 3), 2)) / 3;
            text(x, y, num2str(i), 'FontSize', 10, 'Color', 'black', 'HorizontalAlignment', 'center');
        end

    elseif option == "region"

        % 只显示特定区域附近的编号
        xx = 2;
        yy = 0;
        radius = 0.2;

        for i = 1:ElementsNum

            x = (Nodes(Elements(i, 1), 1) + Nodes(Elements(i, 2), 1) + Nodes(Elements(i, 3), 1)) / 3;
            y = (Nodes(Elements(i, 1), 2) + Nodes(Elements(i, 2), 2) + Nodes(Elements(i, 3), 2)) / 3;

            if abs(x - xx) < radius && abs(y - yy) < radius
                text(x, y, num2str(i), 'FontSize', 10, 'Color', 'black', 'HorizontalAlignment', 'center');
            end

        end

    else

        i = str2double(option);
        x = (Nodes(Elements(i, 1), 1) + Nodes(Elements(i, 2), 1) + Nodes(Elements(i, 3), 1)) / 3;
        y = (Nodes(Elements(i, 1), 2) + Nodes(Elements(i, 2), 2) + Nodes(Elements(i, 3), 2)) / 3;
        text(x, y, num2str(i), 'FontSize', 10, 'Color', 'black', 'HorizontalAlignment', 'center');

    end

end

function NewNodes = GetNewPosition(Nodes, d)
    times = 5e6;
    num = size(Nodes, 1);
    NewNodes = zeros(num, 2);

    for i = 1:num
        NewNodes(i, 1) = Nodes(i, 1) + d(i, 1) * times;
        NewNodes(i, 2) = Nodes(i, 2) + d(i, 2) * times;
    end

    return;
end

function [strain, stress] = GetStrainAndStress(Nodes, Elements, d)

    % 由位移计算应变和应力，应变应力均定义在三角形单元上，而不是节点上
    % strain中按照epsilonX、epsilonY、gammaXY的顺序排列
    % stress中按照sigmaX、sigmaY、tauXY、sigma_1(主应力)、sigma_2的顺序排列
    % 计算时没考虑三角形顺逆时针的影响，因为comsol读出来都是逆时针的所以正确，如果存在顺时针的情况，可能需要修正

    ElementsNum = size(Elements, 1);
    strain = zeros(ElementsNum, 3);
    stress = zeros(ElementsNum, 6);

    % 选项，在GetElementStiffnessMatrix()中也需要修改
    option1 = 2; % 平面应力问题 —— 1，平面应变问题 —— 2
    % 常数
    E = 2E11;
    nu = 0.3;

    if option1 == 1
        % 平面应力问题
        D = E / (1 - nu ^ 2) * [1, nu, 0; nu, 1, 0; 0, 0, (1 - nu) / 2];
    elseif option1 == 2
        E_ = E / (1 - nu ^ 2);
        nu_ = nu / (1 - nu);
        D = E_ / (1 - nu_ ^ 2) * [1, nu_, 0; nu_, 1, 0; 0, 0, (1 - nu_) / 2];
    else
        disp("非法选项")
    end

    for elementIndex = 1:ElementsNum

        xi = Nodes(Elements(elementIndex, 1), 1);
        yi = Nodes(Elements(elementIndex, 1), 2);
        xj = Nodes(Elements(elementIndex, 2), 1);
        yj = Nodes(Elements(elementIndex, 2), 2);
        xk = Nodes(Elements(elementIndex, 3), 1);
        yk = Nodes(Elements(elementIndex, 3), 2);

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

        ui = d(Elements(elementIndex, 1), 1);
        vi = d(Elements(elementIndex, 1), 2);
        uj = d(Elements(elementIndex, 2), 1);
        vj = d(Elements(elementIndex, 2), 2);
        uk = d(Elements(elementIndex, 3), 1);
        vk = d(Elements(elementIndex, 3), 2);

        d_e = [ui; vi; uj; vj; uk; vk];

        epsilon_e = B * d_e;

        sigma_e = D * epsilon_e;

        strain(elementIndex, 1) = epsilon_e(1);
        strain(elementIndex, 2) = epsilon_e(2);
        strain(elementIndex, 3) = epsilon_e(3);

        stress(elementIndex, 1) = sigma_e(1);
        stress(elementIndex, 2) = sigma_e(2);
        stress(elementIndex, 3) = sigma_e(3);

        % 计算主应力
        stress(elementIndex, 4) = (sigma_e(1) + sigma_e(2)) / 2 + sqrt(((sigma_e(1) - sigma_e(2)) / 2) ^ 2 + sigma_e(3) ^ 2);
        stress(elementIndex, 5) = (sigma_e(1) + sigma_e(2)) / 2 - sqrt(((sigma_e(1) - sigma_e(2)) / 2) ^ 2 + sigma_e(3) ^ 2);
        stress(elementIndex, 6) = sqrt(((stress(elementIndex, 4) -stress(elementIndex, 5)) ^ 2 + stress(elementIndex, 4) ^ 2 + stress(elementIndex, 5)) * 1/2);
    end

    return;
end

function [maxStress, belongElement] = GetMaxStress(stress, option)

    ElementsNum = size(stress, 1);
    maxStress = -1;

    if option == "Principal"

        for i = 1:ElementsNum

            if maxStress < stress(i, 4)
                maxStress = stress(i, 4);
                belongElement = i;
            end

        end

    elseif option == "Mises"

        for i = 1:ElementsNum

            mises = stress(i, 6);

            if maxStress < mises
                maxStress = mises;
                belongElement = i;
            end

        end

    end

    return;

end

function F = AddCentralizedLoad(F, nodeIndex, force, option)

    if option == "x"

        F(nodeIndex * 2 - 1) = F(nodeIndex * 2 - 1) + force;

    elseif option == "y"

        F(nodeIndex * 2) = F(nodeIndex * 2) + force;

    else

        disp("无效输入");

    end

end

function OutPutRegionData(Nodes, Elements, xx, yy, radius)

    fprintf('Nodes:\n');
    NodesNum = size(Nodes, 1);

    for i = 1:NodesNum
        x = Nodes(i, 1);
        y = Nodes(i, 2);

        if abs(x - xx) < radius && abs(y - yy) < radius
            fprintf('%d:(%.2f,%.2f)\n', i, x, y);
        end

    end

    fprintf('Elements:\n');
    ElementsNum = size(Elements, 1);

    for i = 1:ElementsNum

        x1 = Nodes(Elements(i, 1), 1);
        x2 = Nodes(Elements(i, 2), 1);
        x3 = Nodes(Elements(i, 3), 1);
        y1 = Nodes(Elements(i, 1), 2);
        y2 = Nodes(Elements(i, 2), 2);
        y3 = Nodes(Elements(i, 3), 2);

        x = (x1 + x2 + x3) / 3;
        y = (y1 + y2 + y3) / 3;

        if abs(x - xx) < radius && abs(y - yy) < radius
            fprintf('%d:(%.2f,%.2f)——(%.2f,%.2f)(%.2f,%.2f)(%.2f,%.2f)\n', i, x, y, x1, y1, x2, y2, x3, y3);
        end

    end

end
