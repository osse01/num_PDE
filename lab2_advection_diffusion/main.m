% Setup
close all
clear
clc

% Declare Variables
a = 1.1;
b = 0.05;
t_end = 0.5;
err = [];
N = [20, 50, 100, 200, 400];
dx_vec = [];

for Nx = N
    dx = 2/(Nx-1);
    % 250 time steps until t_end
    dt = 0.5/250; 
    u_0 = @(x) 1 - sin(pi*x);
    u_l = 1; % alpha
    u_r = 1; % beta
    
    % Help variables
    alpha = a * dt / (4*dx);
    beta  = b * dt / (2*dx^2);
    
    % for the check
    %disp(alpha)
    %rhs = 0.5+beta;
    %disp(rhs)

    % Construct Matrix RHS
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
    % p = plot(x,u_approx);
    % drawnow
    for t = dt:dt:t_end
        u_exact = u_xt(x,t,a,b,50);
        u_approx = updateU(u_approx);
        % p = plot(x,u_approx,x,u_exact);
        % xlim([-1,1]);
        % legend('u aproximate','u exact');
        % title("Wave Plot");
        % drawnow
    end
    e = sqrt(dx)*norm(u_approx - u_exact);
    
    dx_vec = [dx_vec dx];
    err = [err,e];
end

for k = 2:5
    EOC(k) = log(err(k)/err(k-1)) / log(N(k-1)/N(k));
end

% Output For Report
fprintf("\nOutput for Report:\n")
fprintf("--------------------------------\n")
fprintf("Nx \t| dx      \t| E_k\n")
fprintf("--------------------------------\n")
for k = 1:5
        fprintf("%d  \t|%.4f \t|%.4f \n", N(k), dx_vec(k), err(k))
end
fprintf("--------------------------------\n")

plot(N(2:5),EOC(2:5),'-x')
xlim([15,450])
ylim([0,3])
xlabel('N_x');
ylabel('EOC');
title('EOC')
xticks(N);
xticklabels(arrayfun(@num2str, N, 'UniformOutput', false));



%% 2.4
clc

a = 1;
b = 0;
t_end = 2;
Nx = 300;
dx = 3/(Nx-1);
dt = 0.5/250; % 250 time steps until t_end
u_0 = @(x) exp(-1000*(x-0.5).^2);

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

x = linspace(0,3,Nx)';
u_approx = u_0(x);
p = plot(x,u_approx);
drawnow
for t = dt:dt:t_end
    u_approx = updateU(u_approx);
    p = plot(x,u_approx);
    xlim([0,3]);
    legend('u aproximate');
    title("Wave Plot");
    xlabel('x')
    ylabel('u')
    drawnow
end

%% 2.5
a = 1;
t_end = 15;
Nx = 300;
b = 0.01;
dx = 3*b/a;
dt = t_end/250; % 250 time steps until t_end
u_0 = @(x) x;

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

x = linspace(0,1,Nx)';
u_approx = u_0(x);
p = plot(x,u_approx);
drawnow
for t = dt:dt:t_end
    u_approx = updateU(u_approx);
    p = plot(x,u_approx);
    xlim([0,1]);
    legend('u aproximate');
    title("Wave Plot");
    xlabel('x')
    ylabel('u')
    drawnow
end