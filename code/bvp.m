function Z_dot = bvp(t, Z, param)
% DGL Modell des Flugzeugs

    % control
    T = param{1};
    C_L = param{2};
    
    % parameter
    [F, C_D_0, alpha, beta, k, m, g] = param{[3:length(param)]};
    
    % X_dot = f(t, X, ...)
    Z_dot([1:4], 1) = dyn_model(t, Z([1:4, 1]), param);
    
    % lam_dot(t) = -H_x
    Z_dot([5:8], 1) = -[
        0;
        (C_D_0 + k * C_L^2)*Z(3)^2*Z(7)-C_L*(3)*Z(8)*((alpha*beta*F*exp(-beta*Z(2)))/(2*m));
        Z(5)*cos(Z(4))+Z(6)*sin(Z(4))-(C_D_0 + k * C_L^2)*Z(3)^2*Z(7)-C_L*(3)*Z(8)*((alpha*beta*F*exp(-beta*Z(2)))/(2*m))+Z(8)*g*cos(Z(2))/Z(3)^2;
        -Z(5)*sin(Z(4)) + Z(6)*Z(3)*cos(Z(4)) + Z(7)*g*cos(Z(4)) - Z(8)*sin(Z(4));
    ];
    

end
