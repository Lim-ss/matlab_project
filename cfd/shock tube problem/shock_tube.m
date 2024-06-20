function [x, rho, u, p] = shock_tube(nx)
    % 参数设置
    % nx = 1600; % 网格数
    L = 1.0; % 管道长度
    dx = L / nx; % 网格间距
    x = linspace(0, L, nx); % 网格点坐标
    t_end = 0.18; % 终止时间
    dt = 0.0001; % 时间步长
    gamma = 1.4; % 比热比

    % 初始条件
    rho = ones(1, nx) * 0.125;
    u = zeros(1, nx);
    p = ones(1, nx) * 0.1;
    rho(1:round(nx / 2)) = 1.0;
    p(1:round(nx / 2)) = 1.0;

    %用于迭代时临时存放数据
    rho_copy = ones(1, nx);
    u_copy = zeros(1, nx);
    E_copy = zeros(1, nx);

    % 计算初始总能
    E = p / (gamma - 1) + 0.5 * rho .* u .^ 2;

    % 时间推进
    t = 0;

    while t < t_end
        % Roe格式计算通量
        [F1, F2, F3] = roe_flux(rho, u, p, E, gamma, nx);
        [R1, R2, R3] = R(rho, u, p, E, nx, x, dx);

        % 更新变量
        rho_copy(2:end - 1) = rho(2:end - 1) - dt / dx * (F1(2:end - 1) - F1(1:end - 2)) + dt * R1(2:end - 1);
        u_copy(2:end - 1) = (rho(2:end - 1) .* u(2:end - 1) - dt / dx * (F2(2:end - 1) - F2(1:end - 2)) + dt * R2(2:end - 1)) ./ rho_copy(2:end - 1); %注意这个除数需要是新的值
        E_copy(2:end - 1) = E(2:end - 1) - dt / dx * (F3(2:end - 1) - F3(1:end - 2)) + dt * R3(2:end - 1);
        p = (gamma - 1) * (E - 0.5 * rho .* u .^ 2);
        rho(2:end - 1) = rho_copy(2:end - 1);
        u(2:end - 1) = u_copy(2:end - 1);
        E(2:end - 1) = E_copy(2:end - 1);

        % 更新时间
        t = t + dt;
    end

end

function [F1, F2, F3] = roe_flux(rho, u, p, E, gamma, nx)
    % 计算通量
    F1 = rho .* u;
    F2 = rho .* u .^ 2 + p;
    F3 = (E + p) .* u;

    % 计算Roe平均量
    for i = 1:nx - 1
        % 左右状态
        rhoL = rho(i); rhoR = rho(i + 1);
        uL = u(i); uR = u(i + 1);
        pL = p(i); pR = p(i + 1);

        % Roe平均量
        rL = sqrt(rhoL); rR = sqrt(rhoR);
        rhoRoe = rL * rR;
        uRoe = (rL * uL + rR * uR) / (rL + rR);
        hL = (gamma * pL) / ((gamma - 1) * rhoL);
        hR = (gamma * pR) / ((gamma - 1) * rhoR);
        hRoe = (rL * hL + rR * hR) / (rL + rR);
        aRoe = sqrt((gamma - 1) * (hRoe - 0.5 * uRoe ^ 2));

        % Differences
        delta_rho = rhoR - rhoL;
        delta_u = uR - uL;
        delta_p = pR - pL;

        % 特征值
        lambda1 = abs(uRoe);
        lambda2 = abs(uRoe + aRoe);
        lambda3 = abs(uRoe - aRoe);

        % 特征向量
        R1 = [1, uRoe, (uRoe ^ 2) / 2];
        R2 = [1, uRoe + aRoe, hRoe + uRoe * aRoe] * (rhoRoe / 2 / aRoe);
        R3 = [1, uRoe - aRoe, hRoe - uRoe * aRoe] * (-rhoRoe / 2 / aRoe);

        % 幅值
        w1 = delta_rho - delta_p / (aRoe ^ 2);
        w2 = delta_u + delta_p / (rhoRoe * aRoe);
        w3 = delta_u - delta_p / (rhoRoe * aRoe);

        % 通量
        F1_half = 0.5 * (F1(i) + F1(i + 1)) - 0.5 * (lambda1 * w1 * R1(1) + lambda2 * w2 * R2(1) + lambda3 * w3 * R3(1));
        F2_half = 0.5 * (F2(i) + F2(i + 1)) - 0.5 * (lambda1 * w1 * R1(2) + lambda2 * w2 * R2(2) + lambda3 * w3 * R3(2));
        F3_half = 0.5 * (F3(i) + F3(i + 1)) - 0.5 * (lambda1 * w1 * R1(3) + lambda2 * w2 * R2(3) + lambda3 * w3 * R3(3));

        % 更新通量
        F1(i) = F1_half;
        F2(i) = F2_half;
        F3(i) = F3_half;
    end

end

function [R1, R2, R3] = R(rho, u, p, E, nx, x, dx)

    A = 4.2 - 2 .* sqrt(4 - (x - 0.5) .^ 2);
    dA = zeros(1, nx);

    for i = 2:nx - 1
        dA(i) = (A(i) - A(i - 1)) / dx;
    end

    % 计算通量
    R1 = -rho .* u .* dA ./ A;
    R2 =- (rho .* (u .^ 2)) .* dA ./ A;
    R3 =- (E + p) .* u .* dA ./ A;

end
