function [x] = thomasAlg(a,b,c,d)
% Standard implementation of Thomas Algorithm 
%   Solves Ax=d where A is a tridiagonal system with a,b,c diagonals
% Returns the solution x

N = length(b); % size of x
x = zeros(N,1);

for i = 2:N
    w = a(i-1) / b(i-1);
    b(i) = b(i) - w*c(i-1);
    d(i) = d(i) - w*d(i-1);
end
x(N) = d(N) / b(N);
for i=N-1:-1:1
    x(i) = ( d(i) - c(i) * x(i+1) ) / b(i);
end
end