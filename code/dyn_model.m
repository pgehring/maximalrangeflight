function X_dot = dyn_model(t, X, param)
% DGL Modell des Flugzeugs
    
    % control
    T = param{1};
    C_L = param{2};
    
    % parameter
    [F, C_D_0, alpha, beta, k, m, g] = param{[3:length(param)]};

    
    X_dot =         [                                 X(1) * cos(X(4)); % x_dot
                                                      X(1) * sin(X(4)); % h_dot
                        (1/m) * (T - force_model(X(3),X(2), C_L, F, C_D_0, alpha, beta, k, "drag") - m*g*sin(X(4))); % v_dot
                     (1/(m*X(3))) * (force_model(X(3),X(2), C_L, F, C_D_0, alpha, beta, k, "lift") - m*g*cos(X(4)))];% gamma_dot

end
