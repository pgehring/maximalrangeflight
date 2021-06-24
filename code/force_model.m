function F_ret = force_model(v, h, C_L, F, C_D_0, alpha, beta, k, force)
% force function F_ret(v, h, C_L)
    
    rho = alpha * exp(- beta * h);
    q = 0.5 * rho * v^2;
    
    switch force
        case "drag"
            coeff = C_D_0 + k * C_L^2;
        case "lift"
            coeff = C_L;
            
    end

    F_ret = F * coeff * q;

end