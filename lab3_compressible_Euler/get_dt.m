function [dt,lambda_max] = get_dt(u,dx,CFL,gamma)
    prim = cons2prim(u,gamma);
    rho = prim(1,:);
    v = prim(2,:);
    p = prim(3,:);

    a = sqrt(gamma*p./rho);
    lambda_max = max(abs(v) + a);
    dt = CFL*dx/lambda_max;
end