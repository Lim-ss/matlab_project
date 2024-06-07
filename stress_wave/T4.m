% 《应力波基础讲义》
% 塑性中心波的突然卸载

% 绘制X-t图
hold on

plot([-5, 35], [0, 0], "black", 'LineWidth', 1);
plot([0, 0], [0, 20], "black", 'LineWidth', 1);
plot([-4, -4], [0, 3], "black", 'LineWidth', 1);
plot([-4, 0], [3, 3], "black", 'LineWidth', 1);

mplot(0, 0, 30, 0.5, "black");
mplot(0, 0, 20, 1, "red");
mplot(0, 0, 18, 1.05, "red");
mplot(0, 0, 16, 1.1, "red");


A3 = nplot(0, 0, 1.4, 0, 3, 0.5, "red");
A2 = nplot(0, 0, 1.5, 0, 3, 0.5, "red");
A1 = nplot(0, 0, 1.6, 0, 3, 0.5, "red");
mplot(0, 3, A3(1), 0.5, "black");

E1 = mplot(A1(1), A1(2), 0, 0.5, "black");
E2 = mplot(A2(1), A2(2), 0, 0.5, "blue");
E3 = mplot(A3(1), A3(2), 0, 0.5, "green");

B1 = nplot(E1(1), E1(2), 0.5, 0, 0, 1.25, "black");
B2 = nplot(E2(1), E2(2), 0.5, 0, 0, 1.2, "blue");
B3 = nplot(E3(1), E3(2), 0.5, 0, 0, 1.15, "green");
mplot(0, 0, B1(1), 1.25, "red");
mplot(0, 0, B2(1), 1.2, "red");
mplot(0, 0, B3(1), 1.15, "red");

T1 = mplot(B1(1), B1(2), 0, 0.5, "black");
T2 = mplot(B2(1), B2(2), 0, 0.5, "blue");
T3 = mplot(B3(1), B3(2), 0, 0.5, "green");

C1 = nplot(T1(1), T1(2), 0.5, 0, 0, 1.1, "black");
C2 = nplot(T2(1), T2(2), 0.5, 0, 0, 1.05, "blue");
C3 = nplot(T3(1), T3(2), 0.5, 0, 0, 1.0, "green");
mplot(0, 0, C1(1), 1.1, "red");
mplot(0, 0, C2(1), 1.05, "red");
mplot(0, 0, C3(1), 1.0, "red");


function x2y2 = mplot(x1, y1, x2, k, color)

    if x1 < x2
        y2 = y1 + (x2 - x1) * k;
        x2y2(1) = x2;
        x2y2(2) = y2;
        plot([x1, x2], [y1, y2], color, 'LineWidth', 1);

        return;
    else
        y2 = y1 + (x1 - x2) * k;
        x2y2(1) = x2;
        x2y2(2) = y2;
        plot([x1, x2], [y1, y2], color, 'LineWidth', 1);

        return;
    end

end

function x2y2 = nplot(x1, y1, k1, x3, y3, k3, color)

    x2 = ((y3 - y1) + (k1 * x1 - k3 * x3)) / (k1 - k3);
    y2 = y1 + (x2 - x1) * k1;
    plot([x1, x2], [y1, y2], color, 'LineWidth', 1);
    x2y2(1) = x2;
    x2y2(2) = y2;
    return;
end
