function f_star = roe_flux(u_ll, u_rr, gamma)
%
%  Compute the numerical flux from the Roe approximate Riemann solver
%  strategy. This Riemann solver is "sophisticated" as it penalises the
%  three waves (two acoustic and one entropy) with different weightings to
%  ensure upwinding while minimising the amount of artificial dissipation
%
%  INPUT:
%    u_ll: Left solution state vector [density, momentum, energy]
%    u_rr: Right solution state vector [density, momentum, energy]
%    gamma: Specific heat ratio
%
%  OUPUT:
%    f_star: Numerical flux vector using the approximate Riemann solver of Roe
%
    f_star = zeros(3,1);

    % Compute the flux on either side
    f_ll = flux(u_ll, gamma);
    f_rr = flux(u_rr, gamma);

    % get a local copy of the primitive variables
    prim_ll = cons2prim(u_ll, gamma);
    prim_rr = cons2prim(u_rr, gamma);

    % enthalpy
    H_ll = (u_ll(3) + prim_ll(3)) / u_ll(1);
    H_rr = (u_rr(3) + prim_rr(3)) / u_rr(1);

    % compute the Roe averages
    z = sqrt(prim_rr(1) / prim_ll(1));
    rho_hat = z * prim_ll(1);
    v_hat = (prim_ll(2) + z * prim_rr(2)) / (1 + z);
    H_hat = (H_ll + z * H_rr) / (1 + z);

    c2_hat = (gamma-1)*(H_hat - 0.5*v_hat^2);
    c_hat  = sqrt(c2_hat);

    % evaluate right eigenvectors at the Roe states
    r1_hat = [1; v_hat - c_hat; H_hat - v_hat * c_hat];
    r2_hat = [1; v_hat; 0.5 * v_hat^2];
    r3_hat = [1; v_hat + c_hat; H_hat + v_hat * c_hat];

    % evaluate eigenvalues at Roe states
    lambda_hat = [v_hat - c_hat; v_hat; v_hat + c_hat];
    lambda_hat = abs(lambda_hat);

    % compute the jumps from the Roe matrix
    omega_hat = [0.5 * ((prim_rr(3) - prim_ll(3)) - rho_hat*c_hat*(prim_rr(2) - prim_ll(2)));...
                 -(prim_rr(3) - prim_ll(3)) + c2_hat*(prim_rr(1) - prim_ll(1));...
                 0.5 * ((prim_rr(3) - prim_ll(3)) + rho_hat*c_hat*(prim_rr(2) - prim_ll(2)))] ./ c2_hat;

    f_star = 0.5*(f_ll + f_rr) - 0.5*( lambda_hat(1)*omega_hat(1)*r1_hat ...
                                      +lambda_hat(2)*omega_hat(2)*r2_hat ...
                                      +lambda_hat(3)*omega_hat(3)*r3_hat);

end