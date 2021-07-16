function F_ret = force_model(v, h, C_L, param, force)
% force function F_ret(v, h, C_L)
    
    alpha = param.alpha;
    beta = param.beta;
    C_D_0 = param.C_D_0;
    k = param.k;
    F = param.F;
    
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