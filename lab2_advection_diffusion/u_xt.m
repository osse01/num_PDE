function u = u_xt(x, t, a, b, N)
% Exact u as in computer exercise

% prefactor
prefactor = 16 * pi^2 * a * b^3 .* exp((a / (2 * b)) .* (x - (a / 2) * t));

sum1 = 0;
sum2 = 0;
for p = 0:N
    % first
    num1 = (-1)^p * 2 * p .* sin(p * pi * x) .* exp(-p^2 * pi^2 * b * t);
    den1 = a^4 + 8 * (pi * a * b)^2 * (p^2 + 1) + 16 * (pi * b)^4 * (p^2 - 1)^2;
    sum1 = sum1 + num1 ./ den1;

    % second
    k = 2*p + 1;
    num2 = (-1)^p * k .* cos((k/2) * pi * x) .* exp(-(k^2 / 4) * pi^2 * b * t);
    den2 = a^4 + (pi * a * b)^2 * (8*p^2 + 8*p + 10) + ...
           (pi * b)^4 * (4*p^2 + 4*p - 3)^2;
    sum2 = sum2 + num2 ./ den2;
end

u = 1 + prefactor .* (sinh(a / (2 * b)) .* sum1 + cosh(a / (2 * b)) .* sum2);
end
