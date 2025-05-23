function [prim] = cons2prim(u,gamma)
    rho = u(1,:);
    v = u(2,:)./rho;
    p = (gamma-1)*(u(3,:) - rho.*(v.^2)/2);
    prim = [rho;v;p];
end