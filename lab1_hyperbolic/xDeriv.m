function dqdx = xDeriv(q, dx, N)
% Let q be function values evaluated in each point
% Let dx be the resolution and N the number of points
% Lecture => \delta U_j  = (U_{j+1} - U_j{-1})/(2\Delta x)

    dqdx = zeros(size(q)); 
    
    % Interior points 
    for j = 2:N
        dqdx(:,j) = (q(:,j+1) - q(:,j-1)) / (2*dx);
    end
    
    % Left boundary
    dqdx(:,1) = (q(:,2) - q(:,1)) / dx;
    
    % Right boundary
    dqdx(:,N+1) = (q(:,N+1) - q(:,N)) / dx;
end