function X = indirect_model(t, z, params)
% dynamical model of airplane
    [alpha, beta, g, C_D_0, e, F, AR, m, q_max, T, C_L] = params{:};

    k = 1/(pi*e*AR);
    
    X = [   z(4)*sind(z(2));
            1/(2*m*z(4)) * (F*C_L*(z(4))^2*alpha*exp(-beta*z(1)) - 2*m*g*cosd(z(2)));
            z(4)*cosd(z(2));
            1/(2*m) * (2*T + (-C_D_0 - k*C_L^2)*F*z(4)^2*alpha*exp(-beta*z(1)) - 2*m*g*sind(z(2)));
            -(-(alpha* beta* F*exp(- beta*z(1))*C_L*z(4)*z(6))/(2*m) + ((C_D_0+ k*C_L^2)* alpha* beta* F*exp(-beta*z(1))*z(4)^2*z(8))/(2* m));
            -(cosd(z(2))*z(4)*z(5) + ( g*sind(z(2))*z(6))/z(4) - sind(z(2))*z(4)*z(7) - cosd(z(2))* g*z(8));
            0;
            -(sind(z(2))*z(5)+ ((F*alpha*exp(- beta*z(1))*C_L)/(2* m)+( g*cosd(z(2)))/(z(4)^2))*z(6) + cosd(z(2))*z(7) - (( C_D_0+ k*C_L^2)* F* alpha*exp(-beta*z(1))*z(4)*z(8))/ m)];
end