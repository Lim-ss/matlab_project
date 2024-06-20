clear;
clc;

[x,rho, u, p] = shock_tube(1600);

% 绘制结果
figure;
subplot(3, 1, 1);
plot(x, rho);
title('Density');
subplot(3, 1, 2);
plot(x, u);
title('Velocity');
subplot(3, 1, 3);
plot(x, p);
title('Pressure');