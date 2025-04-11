function dqdt = tDeriv(q, dx, t, N, A)
    % q = [u, v]' so u = q[]
    
    % Space derivitives
    dqdx = xDeriv(q, dx, N);
    
    TMP = -1*A*dqdx;
    
    dudt = TMP(1, :);
    dvdt = TMP(2, :);

    dqdt = [dudt; dvdt];


end