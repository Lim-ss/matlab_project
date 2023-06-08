
f = @(a) exp(cos(a));
real_value = integral(f,-pi,pi);
result = zeros(2^10-1);
error = zeros(2^10-1);

for m = 2:2 ^ 10 %子区间数量
    h = 2 * pi / m; %步长
    x = linspace(-pi, pi, m + 1);
    y = f(x);

    result(m - 1) = h / 2 * (f(-pi) + f(pi));

    for i = 1:m - 1
        result(m - 1) = result(m - 1) + h * f(-pi + i * h);
    end

    error(m - 1) = result(m - 1) - real_value;
end
semilogy(2:2 ^ 10, abs(error));
