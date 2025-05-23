function [f] = flux(u,gamma)
    rho = u(1);
    v = u(2);
    p = u(3);
    E = p / (gamma-1) + rho*v^2/2;
    f1 = rho*v;
    f2 = rho*v^2 + p;
    f3 = (E+p)*v;
    f = [f1,f2,f3];
end