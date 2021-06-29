function Z_dot = bvp(t, Z, param)
% DGL Modell des Flugzeugs

    % control
    T = param.T;
    C_L = param.C_L;
    
    % parameter
    F = param.F;
    C_D_0 = param.C_D_0;
    alpha = param.alpha; 
    beta = param.beta;
    k = param.k;
    m = param.m;
    g = param.g;
    
    % X_dot = f(t, X, ...)
    Z_dot([1:4], 1) = dyn_model(t, Z([1:4, 1]), param);
    
    % lam_dot(t) = -H_x
    Z_dot([5:8], 1) = -[
        0;
        Z(7)/(2*m)*(C_D_0 + k * C_L^2)*Z(3)^2*(alpha*beta*F*exp(-beta*Z(2))) - Z(8)/(2*m)*C_L*Z(3)*(alpha*beta*F*exp(-beta*Z(2)));
        Z(5)*cos(Z(4)) + Z(6)*sin(Z(4)) - Z(7)/(m)*(C_D_0 + k * C_L^2)*Z(3)*(alpha*beta*F*exp(-beta*Z(2))) - Z(8)/(2*m)*(C_L*Z(3)*alpha*beta*F*exp(-beta*Z(2))+2*m*g/(Z(3)^2)*cos(Z(2)));
        -Z(5)*Z(3)*sin(Z(4)) + Z(6)*Z(3)*cos(Z(4)) - Z(7)*g*cos(Z(4)) + Z(8)*g/Z(3)*sin(Z(4));
    ];

end
