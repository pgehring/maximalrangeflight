function X_dot = dyn_model(t, X, param)
% DGL Modell des Flugzeugs
    
    % control
    T = param.T;
    C_L = param.C_L;
    
    % parameter
    m = param.m;
    g = param.g;

    
    X_dot =         [                                 X(1) * cos(X(4)); % x_dot
                                                      X(1) * sin(X(4)); % h_dot
                        (1/m) * (T - force_model(X(3),X(2), C_L, param, "drag") - m*g*sin(X(4))); % v_dot
                     (1/(m*X(3))) * (force_model(X(3),X(2), C_L, param, "lift") - m*g*cos(X(4)))];% gamma_dot

end
