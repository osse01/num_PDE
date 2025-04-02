clear 
close all
clc

% Test parameters
N = 100;            % Number of intervals
dx = 1/N;           % Grid spacing
x = linspace(0, 1, N+1)';  % Spatial grid

t1 = ones(5, 1)
t2 = [exp(2); exp(3)]

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

dt1 = xDeriv(t1, dx, N);
disp("Should be 0;")
norm(dt1, 2)

