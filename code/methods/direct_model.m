function X = direct_model(t, z, params)
% dynamical model of airplane
    [alpha, beta, g, C_D_0, e, F, AR, m, q_max, T, C_L] = params{:};

    k = 1/(pi*e*AR);
    
    X = [   z(4)*sind(z(2));
            1/(2*m*z(4)) * (F*C_L*(z(4))^2*alpha*exp(-beta*z(1)) - 2*m*g*cosd(z(2)));
            z(4)*cosd(z(2));
            1/(2*m) * (2*T + (-C_D_0 - k*C_L^2)*F*z(4)^2*alpha*exp(-beta*z(1)) - 2*m*g*sind(z(2)))];
end