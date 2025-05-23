function [cons] = prims2con(prims,gamma)
    u1 = prims(1,:); % rho
    u2 = u1 .* prims(2,:); % rho v
    u3 = prims(3,:)/(gamma-1) + u1.*prims(2,:).^2./2; % E
    cons = [u1;u2;u3];
end