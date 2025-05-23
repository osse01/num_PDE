function [f1,f2,f3] = flux(rho,v,E,gamma)
    p = (gamma - 1)*(E - v^2/2);
    f1 = rho*v;
    f2 = rho*v^2 + p;
    f3 = (E+p)*v;
end