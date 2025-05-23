%% Lab-3
% Setup
close all
clear
clc
display = false;

% Parameters
CFL = 0.95;
gamma = 1.4;
t_end = 0.5;


% Functions
prim_0 = @(x) [2 + 0.5*sin(2*pi*x);
                ones(size(x));
                ones(size(x))]; % Initial data

u_exact = @(x,t) [2 + 0.5*sin(2*pi*(x - t));
                  2 + 0.5*sin(2*pi*(x - t));
                  1 + 0.25*sin(2*pi*(x - t)) + 1/(gamma-1)]; % Exact solution

Nx = [50, 100, 200, 400, 800];

EOC = table('Size', [0 5], ...
            'VariableTypes', {'double', 'double', 'double', 'double', 'double'}, ...
            'VariableNames', {'k', 'Nx', 'dx', 'Ek', 'EOC'});
k = 1;

for J = Nx
    if(display)
        figure(J)
    end
    dx = 1/J;
    x = (0.5 + (0:J-1))*dx;
    u = prims2con(prim_0(x),gamma);
    t = 0;
    while t < t_end
        [dt,lambda_max] = get_dt(u,dx,CFL,gamma);
        t = t + dt; % Time Update
        
        if(true)
            % Rusanov Flux Method
            u_new(:,1) = u(:,1) - (dt/dx) * ...
                (rusanov_flux(u(:,1),u(:,2),gamma,lambda_max) - ...
                rusanov_flux(u(:,J),u(:,1),gamma,lambda_max));
    
            for j = 2:J-1
                u_new(:,j) = u(:,j) - dt/dx * ...
                (rusanov_flux(u(:,j),u(:,j+1),gamma,lambda_max) - ...
                rusanov_flux(u(:,j-1),u(:,j),gamma,lambda_max));
            end
            
            u_new(:,J) = u(:,J) - dt/dx * ...
                (rusanov_flux(u(:,J),u(:,1),gamma,lambda_max) - ...
                rusanov_flux(u(:,J-1),u(:,J),gamma,lambda_max));
        
            u = u_new;
        else
            % Roe Flux Method
            u_new(:,1) = u(:,1) - dt/dx * ...
                (roe_flux(u(:,1),u(:,2),gamma) - ...
                roe_flux(u(:,J),u(:,1),gamma));
            
            for j = 2:J-1
                u_new(:,j) = u(:,j) - dt/dx * ...
                (roe_flux(u(:,j),u(:,j+1),gamma) - ...
                roe_flux(u(:,j-1),u(:,j),gamma));
            end
            
            u_new(:,J) = u(:,J) - dt/dx * ...
                (roe_flux(u(:,J),u(:,1),gamma) - ...
                roe_flux(u(:,J-1),u(:,J),gamma));

            u = u_new;
        end

        u_xt = u_exact(x,t);
        
        % Plot Solution
        if (display)
            plot(x,u(1,:),x,u_xt(1,:));
            legend('u_approx','u_exact')
            drawnow;
        end
        
    end

    % Calculate Error
        u_xt = u_exact(x, t);
        err = u(1,:) - u_xt(1,:);
        Ek = sqrt(sum(err.^2) * dx);
    
        if k > 1
            EOC_k = log(Ek_prev / Ek) / log(dx_prev / dx);
        else
            EOC_k = NaN;
        end
    
        newRow = {k, J, dx, Ek, EOC_k};
        EOC = [EOC; newRow];
    
        % Update previous error and dx
        Ek_prev = Ek;
        dx_prev = dx;
        k = k + 1;
end

%% Exercise 2.3
close all
clear
clc
display = true;
gamma = 1.4;
CFL = 0.95;
t_end = 0.2;
data = readmatrix('sod_reference_solution.txt', 'FileType', 'text', 'CommentStyle', '%');
x_exact = data(:,1);
u_exact = data(:,2:4);
k = 1;

J = 1000;
x = linspace(0,1,J);
u = prims2con(initial(x),gamma);
dx = 1/J;
t = 0;
while t < t_end
    [dt,lambda_max] = get_dt(u,dx,CFL,gamma);
    t = t + dt; % Time Update
    
    if(true)
        % Rusanov Flux Method
        u_new(:,1) = u(:,1) - (dt/dx) * ...
            (rusanov_flux(u(:,1),u(:,2),gamma,lambda_max) - ...
            rusanov_flux(u(:,1),u(:,1),gamma,lambda_max));

        for j = 2:J-1
            u_new(:,j) = u(:,j) - dt/dx * ...
            (rusanov_flux(u(:,j),u(:,j+1),gamma,lambda_max) - ...
            rusanov_flux(u(:,j-1),u(:,j),gamma,lambda_max));
        end
        
        u_new(:,J) = u(:,J) - dt/dx * ...
            (rusanov_flux(u(:,J),u(:,J),gamma,lambda_max) - ...
            rusanov_flux(u(:,J-1),u(:,J),gamma,lambda_max));
        
        u = u_new;
    else
        % Roe Flux Method
        u_new(:,1) = u(:,1) - dt/dx * ...
            (roe_flux(u(:,1),u(:,2),gamma) - ...
            roe_flux(u(:,J),u(:,1),gamma));
        
        for j = 2:J-1
            u_new(:,j) = u(:,j) - dt/dx * ...
            (roe_flux(u(:,j),u(:,j+1),gamma) - ...
            roe_flux(u(:,j-1),u(:,j),gamma));
        end
        
        u_new(:,J) = u(:,J) - dt/dx * ...
            (roe_flux(u(:,J),u(:,1),gamma) - ...
            roe_flux(u(:,J-1),u(:,J),gamma));

        u = u_new;
    end
    
    % Plot Solution
    if (display)
        plot(x,u(1,:));
        hold on;
        plot(x_exact,u_exact(:,1));
        hold off;
        legend('u_approx','u_exact')
        drawnow;
    end
    
end


function prim_0 = initial(x)
N = numel(x);
rho = zeros(1,N);
p = zeros(1,N);
v = zeros(1,N);
    for idx = 1:N
        if (x(idx)<0.5)
            rho(idx) = 1;
            p(idx) = 1;
        else
            rho(idx) = 0.125;
            p(idx) = 0.1;
        end
    end
prim_0 = [rho;v;p];
end