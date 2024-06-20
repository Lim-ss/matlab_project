clear;
clc;

[x_ref, rho_ref, u_ref, p_ref] = shock_tube(1600);

d_log = zeros(1, 1);
error_log = zeros(1, 1);
m_error=0;

for i = 1:8
    n = 100 * i;
    [x, rho, u, p] = shock_tube(n);

    times = 1600 / n;

    rho_ = zeros(n, 1);
    u_ = zeros(n, 1);
    p_ = zeros(n, 1);

    for j = 1:n
        rho_(j) = rho_ref(floor(times * j));
        u_(j) = u_ref(floor(times * j));
        p_(j) = p_ref(floor(times * j));
    end

    %计算误差积分
    for j = 1:n
        m_error = m_error + (rho(j) - rho_(j)) ^ 2;
    end

    m_error = m_error / n;

    d_log(i) = log10(1/n);
    error_log(i) = log10(m_error);

end

plot(d_log, error_log);
xlabel('log_{10}(d)');
ylabel('log_{10}(error)');
title('误差随网格间距变化');

% % 指定Excel文件路径和文件名
% filename = 'D:/QQ消息记录等数据/757001674/FileRecv/参考解.xlsx';
% rho_ref = readmatrix(filename, 'Range', 'B2:B201');
% u_ref = readmatrix(filename, 'Range', 'C2:C201');
% p_ref = readmatrix(filename, 'Range', 'D2:D201');
