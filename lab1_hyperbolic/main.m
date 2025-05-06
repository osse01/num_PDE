clear 
close all
clc

%% Excercise 2.1
disp("Excercise 2.1")

% Params
N = 100;            
dx = 1/N;           
x = linspace(0, 1, N+1);  

% Uniform partition of [0,1] wih N partitions
f = @(x) exp(2*x);
q = [f(x);zeros(size(x))];


% Demonstrate that function works.
dq = xDeriv(q, dx, N);

figure;
plot(x, dq, 'r')  
hold on
plot(x, q, 'g')  
plot(x, (dq-q)./q, 'b')
xlabel('x');
ylabel('Values');
legend('Derivitive of exp(2x)', 'exp(2x)', 'Difference divided by e^2x (should be 1)', Location='best');
title('Controll of dqdx function');
grid on;

%% Excercise 2.2
disp("Excercise 2.2")

% Params
N = 100;
dx = 1/N;
x = linspace(0,1,N+1);
t = 0;


u = sin(2*pi*x);
v = cos(2*pi*x);
q = nan(2,(N+1));
q(1,:) = u;  
q(2,:) = v;

A = [1, 2; 2, 1];

TMP = tDeriv(q, dx, t, N, A);
dudt = TMP(1, :);
dvdt = TMP(2, :);

figure
subplot(2,1,1)
plot(x, dvdt, 'b', x, dudt, 'r');
legend('dv/dt = -2\pi cos(2\pi x)', 'du/dt = 2\picos(2\pi x)', Location='best');
title('dv/dt and du/dt');

subplot(2,1,2)
plot(x, v, 'b', x, u, 'r');
legend('v = cos(2\pi x) ', 'u = sin(2\pi x)', Location='best');
title('v and u');

%% Exercise 2.3
clc
disp("Excercise 2.3")


N = 100;
T = 100;
dt = 1/T;
dx = 1/N;          
x = linspace(0, 1, N+1);
f = @(x) exp(-100*(x-0.4).^2);
q = [f(x);zeros(size(x))];

figure(4)
p = plot(q(1,:));

for t = 1:T
    q = stepByRK3(q,t,dx,dt,N,A);
    p = plot(q(1,:));
    title("Wave Plot");
    drawnow;
end

%% Exercise 2.4/5
clc
disp("Excercise 2.4")

A = [1, 2; 2, 1];
N = 100;
dx = 1/N;
lambda_max = max(eig(A));
CFL = 1;
dt = CFL*dx/lambda_max; % \Delta t = \frac{\text{CFL} \Delta x}{|\lambda_{\text{max}}}


q_exact = @(x, t) 0.5 * [exp(-100*(x - 3*t - 0.4).^2) + exp(-100*(x + t - 0.4).^2) - exp(-100*(x/3 - t + 0.4).^2); 
                   exp(-100*(x - 3*t - 0.4).^2) - exp(-100*(x + t - 0.4).^2) - exp(-100*(x/3 - t + 0.4).^2)];
x = linspace(0, 1, N+1);
f = @(x) exp(-100*(x-0.4).^2);
q = [f(x);
    zeros(size(x))];

t = 0;
t_end = 0.75;
figure(5)
[P,E] = eig(A);
Q_plot = [];
timeVec = [];
while t < t_end
    if t + dt > t_end
        dt = t_end - t;
    end
    t = t + dt;
    q = stepByRK3(q,t,dx,dt,N,A);
    
    % q(1,1) = 0;
    % w = P \ q(:,end);  % Characteristic variables at x = 0
    % w(1) = 0; % w^- (1,t)
    % q(:,end) = P * w;

    q_e = q_exact(x,t);
    Q_plot = [Q_plot; q(1,:)];
    timeVec = [timeVec, t];
    p = plot(x,q(1,:),'r',x,q(2,:),'r',x,q_e(1,:),'b', x,q_e(2,:),'b');
    ylim([-0.5,1]);
    legend('q aproximate','q exact');
    title("Wave Plot");
    drawnow;
end
surf(x,timeVec,Q_plot)
xlabel('x')
ylabel('t')

%% Exercise 2.6
clc
A = [1, 2; 2, 1];

lambda_max = max(eig(A));
[P,E] = eig(A);
CFL = 1;


q_exact = @(x, t) 0.5 * [exp(-100*(x - 3*t - 0.4).^2) + exp(-100*(x + t - 0.4).^2) - exp(-100*(x/3 - t + 0.4).^2); 
                   exp(-100*(x - 3*t - 0.4).^2) - exp(-100*(x + t - 0.4).^2) - exp(-100*(x/3 - t + 0.4).^2)];
f = @(x) exp(-100*(x-0.4).^2);
Error = [];
N_values = [50, 100, 200, 400, 800];
for N=N_values
    x = linspace(0, 1, N+1);
    q = [f(x);
        zeros(size(x))];
    dx = 1/N;
    dt = CFL*dx/lambda_max; % \Delta t = \frac{\text{CFL} \Delta x}{|\lambda_{\text{max}}}
    
    t = 0;
    t_end = 0.75;
    e = 0;
    while t < t_end
        if t + dt > t_end
            dt = t_end - t;
        end
        t = t + dt;
        q = stepByRK3(q,t,dx,dt,N,A);
        % q(1,1) = 0;
        % w = P \ q(:,end);  % Characteristic variables at x = 0
        % w(1) = 0; % w_minus(1,t)
        % q(:,end) = P * w;

        q_e = q_exact(x,t);
    end
    e = q - q_e;
    E_j = sqrt(e(1,:).^2 + e(2,:).^2);
    error = sqrt(dx * sum(E_j.^2));
    Error = [Error; error];
end
loglog([50, 100, 200, 400, 800], Error,'o-')
hold on;
p = polyfit(log(N_values), log(Error), 1);
fitted_line = exp(polyval(p, log(N_values)));
loglog(N_values, fitted_line, '--','DisplayName', sprintf('Fit (p=%.2f)', p(1)));
legend('Error',sprintf('Fit (p=%.2f)', p(1)))
xlabel('N');
ylabel('Error');
title('Error Convergence: RK3 Method for Linear System (CFL = 1)')
xticks(N_values);
xticklabels(arrayfun(@num2str, N_values, 'UniformOutput', false));