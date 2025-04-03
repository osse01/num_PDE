function dqdt = tDeriv(q, dx, t, N, A)
    % q = [u, v]' so u = q[]
    
    % Space derivitives
    dudx = xDeriv(q(:, 1), dx, N);
    dvdx = xDeriv(q(:, 2), dx, N);
    
    TMP = -1*A*[dudx, dvdx]';
    
    dudt = TMP(1, :);
    dvdt = TMP(2, :);


    dqdt = [dudt; dvdt];


end