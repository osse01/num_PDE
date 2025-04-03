clear 
close all
clc

%% Excercise 2.1
% Params
N = 100;            
dx = 1/N;           
x = linspace(0, 1, N+1)';  

% Uniform partition of [0,1] wih N partitions
f = @(x) exp(2*x);
q = f(x);

function dqdx = xDeriv(q, dx, N)
% Let q be function values evaluated in each point
% Let dx be the resolution and N the number of points
% Lecture => \delta U_j  = (U_{j+1} - U_j{-1})/(2\Delta x)

    dqdx = zeros(size(q)); 
    
    % Interior points 
    for j = 2:N
        dqdx(j) = (q(j+1) - q(j-1)) / (2*dx);
    end
    
    % Left boundary
    dqdx(1) = (q(2) - q(1)) / dx;
    
    % Right boundary
    dqdx(N+1) = (q(N+1) - q(N)) / dx;
end
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