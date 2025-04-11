function [q] = stepByRK3(q,t,dx,dt,N,A)
% q = [u,v]'
% t evaluation time
% dx space step
% dt time step
% N Resolution
a = [0, -5/9, -153/128];
b = [0, 1/3, 3/4];
g = [1/3, 15/16, 8/15];

G = zeros(size(q));
for k = 1:3
    tk = t + dt*b(k);
    dqdt = tDeriv(q,dx,tk,N,A);
    G = a(k)*G + dqdt;
    q = q + g(k)*dt*G;
end