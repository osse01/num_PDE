clear 
close all
clc

%% Excercise 2.1
disp("Excercise 2.1")

% Params
N = 100;            
dx = 1/N;           
x = linspace(0, 1, N+1)';  

% Uniform partition of [0,1] wih N partitions
f = @(x) exp(2*x);
q = f(x);


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
x = linspace(0,1,N+1)';
t = 0;


u = sin(2*pi*x);
v = cos(2*pi*x);
q = nan((N+1),2);
q(:, 1) = u;  
q(:, 2) = v;

A = [1, 0; 0, 1];

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