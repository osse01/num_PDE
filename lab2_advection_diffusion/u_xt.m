function u = u_xt(x, t, a, b, N)
% u_xt Evaluates the series solution u(x,t) for given constants a, b
%
%   u = u_xt(x, t, a, b, N)
%   Computes the solution for spatial coordinate x, time t, constants a and b,
%   and the number of series terms N (integer >= 0).
%
%   Inputs:
%       x - spatial coordinate (scalar or array)
%       t - time (scalar or array)
%       a - constant (scalar)
%       b - constant (scalar)
%       N - number of terms in the series (integer)
%
%   Output:
%       u - solution u(x,t)

% Precompute the exponential prefactor
prefactor = 1 + 16 * pi^2 * a * b^3 .* exp((a / (2 * b)) .* (x - (a / 2) * t));

% Compute the two series components
sum1 = 0;
sum2 = 0;

for p = 0:N
    % ----- First Series -----
    num1 = (-1)^p * 2 * p .* sin(p * pi * x) .* exp(-p^2 * pi^2 * b * t);
    den1 = a^4 + 8 * (pi * a * b)^2 * (p^2 + 1) + 16 * (pi * b)^4 * (p^2 - 1)^2;
    sum1 = sum1 + num1 ./ den1;

    % ----- Second Series -----
    k = 2*p + 1;
    num2 = (-1)^p * k .* cos((k/2) * pi * x) .* exp(-(k^2 / 4) * pi^2 * b * t);
    den2 = a^4 + (pi * a * b)^2 * (8*p^2 + 8*p + 10) + ...
           (pi * b)^4 * (4*p^2 + 4*p - 3)^2;
    sum2 = sum2 + num2 ./ den2;
end

% Final assembly of the full solution
u = prefactor .* (sinh(a / (2 * b)) .* sum1 + cosh(a / (2 * b)) .* sum2);
end
