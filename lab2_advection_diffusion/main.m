%% Main

% Init
a = 1.1;
b = 0.05;
t_end = 0.5;
Nx = 50; % 20, 50, 100, 200, 400.
dx = 2/(Nx-1);
dt = 0.5/250; % 250 time steps until t_end
u_0 = @(x) 1 - sin(pi*x);
u_l = 1; % alpha
u_r = 1; % beta

% Help variables
alpha = a * dt / (4*dx);
beta  = b * dt / (2*dx^2);

% Construct matrix RHS
main_diag = ones(Nx, 1) * (1 - 2 * beta);
lower_diag = ones(Nx-1, 1) * (alpha + beta);
upper_diag = -ones(Nx-1, 1) * (alpha - beta);
main_diag(1) = 1;
upper_diag(1) = 0;
lower_diag(end) = 0;
main_diag(end) = 1;
RHS = diag(main_diag) + ...
      diag(lower_diag, -1) + ...
      diag(upper_diag, 1);

updateU = @(u) thomasAlg([-(alpha + beta)*ones(1,Nx-2),0], ...
                        [1,(1+2*beta)*ones(1,Nx-2),1], ...
                        [0,(alpha-beta)*ones(1,Nx-2)], ...
                        RHS*u);

x = linspace(-1,1,Nx)';
u_approx = u_0(x);
p = plot(x,u_approx);
xlim([-1,1]);
legend('u aproximate','u exact');
title("Wave Plot");
drawnow
for t = 0:dt:t_end
    u_exact = u_xt(x,t,a,b,25);
    u_approx = updateU(u_approx);
    p = plot(x,u_approx,x,u_exact);
    xlim([-1,1]);
    legend('u aproximate','u exact');
    title("Wave Plot");
    drawnow
    
end
