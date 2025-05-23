function [f] = flux(u,gamma)
% u is in conservative variables
    prim = cons2prim(u,gamma);
    
    rho = prim(1,:);
    v = prim(2,:);
    p = prim(3,:);
    E = p/(gamma-1) + 0.5*rho.*v.^2;

    f1 = rho.*v;
    f2 = rho.*(v.^2) + p;
    f3 = (E+p).*v;
    
    % f1 = u(2,:);
    % f2 = (3-gamma)*(u(2,:).^2)./(2*u(1,:)) + (gamma-1)*u(3,:);
    % f3 = gamma*u(3,:).*u(2,:)./u(1,:) - (gamma-1)*(u(2,:).^3)./(2*u(1,:).^2);
    f = [f1;f2;f3];
end