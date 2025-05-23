function [F_star] = rusanov_flux(u_ll,u_rr,gamma,lambda_max)
%  INPUT:
%    u_ll: Left solution state vector [density, momentum, energy]
%    u_rr: Right solution state vector [density, momentum, energy]
%    gamma: Specific heat ratio
%
%  OUPUT:
%    f_star: Numerical flux vector using the approximate Riemann solver of
%    Rusanov
% 
% lambda_max is a local estimate of the maximum wave speed
    
F_star = ( flux(u_rr,gamma) + flux(u_ll,gamma) )/2 - ...
        lambda_max * (u_rr-u_ll) / 2;
end