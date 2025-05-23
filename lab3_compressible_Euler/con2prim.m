function [u1,u2,u3] = con2prim(rho, v, E)
    u1 = rho;
    u2 = rho*v;
    u3 = E;
end