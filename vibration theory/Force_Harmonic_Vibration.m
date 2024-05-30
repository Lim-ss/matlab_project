%《振动理论及应用》（第五版）P82 题3.6
% 绘制受迫振动响应的图像

zeta = 0.01; % 阻尼系数
ratio = 0:0.1:5; % omega/omega_n

realPart = (1 - (ratio) .^ 2) ./ ((1 - (ratio) .^ 2) .^ 2 + (2 * zeta * ratio) .^ 2);
imaginaryPart =- (2 * zeta * ratio) ./ ((1 - (ratio) .^ 2) .^ 2 + (2 * zeta * ratio) .^ 2);

% 第一个子图
subplot(1, 2, 1);
plot(ratio, realPart, 'r');
hold on;
plot(ratio, imaginaryPart, 'g');
title('zeta = 0.01');
xlabel('w/wn');
ylabel('H');
xlim([0, 5]); % 限制显示范围
ylim([-10, 10]);
hold off;

zeta = 0.02;
ratio = 0:0.01:5;

realPart = (1 - (ratio) .^ 2) ./ ((1 - (ratio) .^ 2) .^ 2 + (2 * zeta * ratio) .^ 2);
imaginaryPart =- (2 * zeta * ratio) ./ ((1 - (ratio) .^ 2) .^ 2 + (2 * zeta * ratio) .^ 2);

% 第二个子图
subplot(1, 2, 2);
plot(ratio, realPart, 'r');
hold on;
plot(ratio, imaginaryPart, 'g');
title('zeta = 0.02');
xlabel('w/wn');
ylabel('H');
xlim([0, 5]); % 限制显示范围
ylim([-10, 10]);
hold off;
