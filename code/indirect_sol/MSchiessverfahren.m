%% Multiple shooting method
function [t_out,y_out,i,eta,Norm_F] = MSchiessverfahren(tspan,eta_0,g,g_y,r,r_y_a,r_y_b,ode_method,AbsTol,RelTol,StopTol,maxit)
    % Step size control
    beta = 0.5;
    sigma = 0.01;
    % ODE Parameters
    options = odeset('RelTol',RelTol,'AbsTol',AbsTol);

    N = size(eta_0,2);

    % Parameter und Einstellungen
    eta = eta_0;
    % Zeitinterval
    t = linspace(tspan(1),tspan(2),N+1)';
    Norm_F = [];

    % Problemgröße
    n_y = size(eta_0,1);
    n_r = length(r(eta_0,eta_0));
    % Allokate Memory
    F = zeros((N-1)*n_y+n_r,1);
    F_jac = zeros((N-1)*n_y+n_r,N*n_y);
    S_a = eye(n_y); % Einheitsmatrix als Startwert für Sensitivitäts-DGL
    S_b_AWP = zeros(n_y);

    function F_Fjac_func(eta)
        % Zurücksetzen der alten Werte
        t_out = [];
        y_out = [];
        % AWP lösen und aufstellen von F(eta) und der Jacobimatrix
        for j = 1:N
            % ODE und Sensitivitäts DGL lösen 
            sol_y = ode_method(g,[t(j),t(j+1)],eta(:,j),options);
            y_AWP = sol_y.y';
            t_AWP = sol_y.x';
            S_t =@(t,S) reshape(g_y(t,deval(sol_y,t))*reshape(S,n_y,n_y),[],1); % Marix-Matrix-Produnkt und zurück in Spaltenvektor 
            sol_S = ode_method(S_t,[t(j),t(j+1)],reshape(S_a,[],1),options);
            S_b_AWP(:,:) = reshape(sol_S.y(:,end)',n_y,n_y); % Spaltenvektor zu Matrix
            % Bilden von F(eta) und der Jacobimatrix
            if j ~= N
                F_jac(((j*n_y)-(n_y-1)):(j*n_y),((j*n_y)-(n_y-1)):(j*n_y)) = S_b_AWP;
                F(((j*n_y)-(n_y-1)):(j*n_y),1)= y_AWP(end,:)' - eta(:,j+1);
                F_jac(((j*n_y)-(n_y-1)):(j*n_y),(((j+1)*n_y)-(n_y-1)):((j+1)*n_y)) = - eye(n_y);
            else
                F((end-n_r+1):end,:) = r(eta(:,1),y_AWP(end,:));
                F_jac(((j*n_y)-(n_y-1)):((j-1)*n_y)+n_r,((j*n_y)-(n_y-1)):(j*n_y)) = r_y_b(eta(:,1),y_AWP(end,:)) * S_b_AWP;
                F_jac((end-n_r+1):end,1:n_y) = r_y_a(eta(:,1),y_AWP(end,:));
            end
            % Ausgabe abspeichern
            t_out((end+1):(end+size(t_AWP,1)),:) = t_AWP;
            y_out((end+1):(end+size(t_AWP,1)),:) = y_AWP;
        end
    end

    % Start des Mehrschiesverfahrens
    i = 0;
    while (i < maxit)
        % AWP lösen und aufstellen von F(eta) und der Jacobimatrix
        F_Fjac_func(eta);
        % Abbruchbedingung
        Norm_F(i+1) = norm(F);
        if (norm(F) <= StopTol)
            break; 
        end
        % Berechnung der Newton-Richtung d 
        d = - F_jac \ F;
        d_r = reshape(d,n_y,N);
        % Armijo-Regel zur Schrittweitensteuerung
        t_armijo = 1; % initial t
        F_Fjac_func(eta + t_armijo * d_r);
        while (F > (F + sigma * t_armijo * F_jac * d))
            t_armijo = beta*t_armijo;
            F_Fjac_func(eta + t_armijo * d_r);
        end
        % Neue Werte eta und i
        eta = eta + t_armijo*d_r;
        i = i + 1;
    end
end



% %% Multiple shooting method
% function [t_out,y_out,i,eta,Norm_F] = MSchiessverfahren(tspan,eta_0,g,g_y,r,r_y_a,r_y_b,ode_method,AbsTol,RelTol,StopTol,maxit)
%     % Step size control
%     beta = 0.5;
%     sigma = 0.01;
%     % ODE Parameters
%     options = odeset('RelTol',RelTol,'AbsTol',AbsTol);
% 
%     N = size(eta_0,2);
% 
%     % Parameter und Einstellungen
%     eta = eta_0;
%     % Zeitinterval
%     t = linspace(tspan(1),tspan(2),N+1)';
%     Norm_F = [];
% 
%     % Problemgröße
%     n_y = size(eta_0,1);
%     n_r = length(r(eta_0,eta_0));
%     % Allokate Memory
%     eta_before = zeros(N*n_y,1);
%     F = zeros((N-1)*n_y+n_r,1);
%     F_before = zeros((N-1)*n_y+n_r,1);
%     F_jac = zeros((N-1)*n_y+n_r,N*n_y);
%     F_jac_before = zeros((N-1)*n_y+n_r,N*n_y);
%     S_a = eye(n_y); % Einheitsmatrix als Startwert für Sensitivitäts-DGL
%     S_b_AWP = zeros(n_y);
% 
%     function F_Fjac_func(eta)
%         % Zurücksetzen der alten Werte
%         t_out = [];
%         y_out = [];
%         % AWP lösen und aufstellen von F(eta) und der Jacobimatrix
%         for j = 1:N
%             % ODE und Sensitivitäts DGL lösen 
%             sol_y = ode_method(g,[t(j),t(j+1)],eta(:,j),options);
%             y_AWP = sol_y.y';
%             t_AWP = sol_y.x';
%             % Bilden von F(eta) und der Jacobimatrix
%             if j ~= N
%                 F(((j*n_y)-(n_y-1)):(j*n_y),1)= y_AWP(end,:)' - eta(:,j+1);
%             else
%                 F((end-n_r+1):end,:) = r(eta(:,1),y_AWP(end,:));
%             end
%             % Ausgabe abspeichern
%             t_out((end+1):(end+size(t_AWP,1)),:) = t_AWP;
%             y_out((end+1):(end+size(t_AWP,1)),:) = y_AWP;
%         end
%         
%         d_b = eta(:) - eta_before;
%         z_b = F - F_before;
%         F_jac = F_jac_before + ((z_b - F_jac_before*d_b)*d_b')/(d_b'*d_b);
%         
%         eta_before = eta(:);
%         F_jac_before = F_jac;
%         F_before = F;
%     end
% 
%     % Start des Mehrschiesverfahrens
%     i = 0;
%     while (i < maxit)
%         % AWP lösen und aufstellen von F(eta) und der Jacobimatrix
%         F_Fjac_func(eta);
%         % Abbruchbedingung
%         Norm_F(i+1) = norm(F);
%         if (norm(F) <= StopTol)
%             break; 
%         end
%         % Berechnung der Newton-Richtung d 
%         d = - F_jac \ F;
%         d_r = reshape(d,n_y,N);
%         % Armijo-Regel zur Schrittweitensteuerung
%         t_armijo = 1; % initial t
%         F_Fjac_func(eta + t_armijo * d_r);
%         while (F > (F + sigma * t_armijo * F_jac * d))
%             t_armijo = beta*t_armijo;
%             F_Fjac_func(eta + t_armijo * d_r);
%         end
%         % Neue Werte eta und i
%         eta = eta + t_armijo*d_r;
%         i = i + 1;
%     end
% end