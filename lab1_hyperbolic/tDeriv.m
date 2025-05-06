function dqdt = tDeriv(q, dx, t, N, A)
    % q = [u, v]' so u = q[]
    
    % Space derivitives
    dqdx = xDeriv(q, dx, N);
    
    TMP = -1*A*dqdx;
    
    dudt = TMP(1, :);
    dvdt = TMP(2, :);

    dudt(1) = 0;
    dvdt(1) = 3*dqdx(2,1);
    
    dvdt(end) = -3*dqdx(2,end);
    dudt(end) = -3*dqdx(1,end);

    dqdt = [dudt; dvdt];


end