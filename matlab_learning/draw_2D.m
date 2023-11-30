x = 0:0.01:2*pi;
y = sin(x);
figure%建立幕布
plot(x,y)
title("y = sin(x)");
xlabel("x");
ylabel("sin(x)");
xlim([0 2*pi]);%视图坐标限制

x = 0:0.01:20;
y1 = 200 * exp(-0.05 * x) .* sin(x);
y2 = -0.8 * exp(-0.5 * x) .* sin(10 * x);
figure;
plot(x,y1);
plot(x,y2);
