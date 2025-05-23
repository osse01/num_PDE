%% Lab-3
CFL = 0.95;
gamma = 1.4;

prim_0 = @(x) [2 + 0.5*sin(2*pi*x);
                ones(size(x));
                ones(size(x))]; % Initial data
t_end = 0.5;
u_exact = @(x,t) [2 + 0.5*sin(2*pi*(x - t));
                  2 + 0.5*sin(2*pi*(x - t));
                  1 + 0.25*sin(2*pi*(x - t)) + 1/(gamma-1)];
Nx = [50, 100, 200, 400, 800];
for J = Nx
    figure(J)
    dx = 1/J;
    x = (0.5 + 0:J-1)*dx;
    u = prims2con(prim_0(x),gamma);
    t = 0;
    while t < t_end
        [dt,lambda_max] = get_dt(u,dx,CFL,gamma);
        for j = 2:J-1
            u(:,j) = u(:,j) - dt/dx * (rusanov_flux(u(:,j),u(:,j+1),gamma,lambda_max) - ...
                                    rusanov_flux(u(:,j-1),u(:,j),gamma,lambda_max));
        end
        u(:,1) = u(:,1) - dt/dx * (rusanov_flux(u(:,1),u(:,2),gamma,lambda_max) - ...
                                    rusanov_flux(u(:,J),u(:,1),gamma,lambda_max));
        u(:,J) = u(:,J) - dt/dx * (rusanov_flux(u(:,J),u(:,1),gamma,lambda_max) - ...
                                    rusanov_flux(u(:,J-1),u(:,J),gamma,lambda_max));
        plot(x,u(1,:));
        hold on
        plot(x,u(2,:));
        plot(x,u(3,:));
        hold off
        legend('u_1','u_2','u_3')
        drawnow;
        t = t + dt;
    end
end
