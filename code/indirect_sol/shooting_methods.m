% shooting_methods.m:
% Klasse an verschiedenen Einschrittverfahren zur Lösung von ODE's

% Date:         27.08.2021
% Author:       Gehring, Philipp / Karus, Heiko / Götz, Felix

classdef shooting_methods
    properties
        %% Variablen für die Berechnung
        ode_method;
        % Step size control
        armijo;
        beta;
        sigma;
        % ODE Parameters
        options;
        % Parameters
        maxit;
        h_min;
        flag;
        StopTol;
        StopTolArmijo;
    end
    methods
        %% Constructor
        function obj = shooting_methods(ode_method,h_min,AbsTol,RelTol,StopTol,StopTolArmijo,maxit,flag)
            obj.ode_method = ode_method;
            % Step size control
            obj.armijo = true;
            obj.beta = 0.5;
            obj.sigma = 0.01;
            % ODE Parameters
            obj.options = odeset('RelTol',RelTol,'AbsTol',AbsTol);
            % Parameters
            obj.maxit = maxit;
            obj.h_min = h_min;
            obj.flag = flag;
            obj.StopTol = StopTol;
            obj.StopTolArmijo = StopTolArmijo;
        end
        
        %% Single shooting method
        function [t_out,y_out,eta,i,Norm_F] = Einfachschiessverfahren(obj,tspan,eta_0,g,g_y,r,r_y_a,r_y_b)
            eta = eta_0(:);
            n_y = length(eta_0);
            n_r = length(r(eta,eta));
            S_b = [];
            %
            Norm_F = [];
            %
            i = 0;
            switch obj.flag
                case 'SensDGL'
                    %% Benötigten Ableitungen mit Sensitivitätsdifferentialgleichung lösen
                    S_a = eye(n_y); % Einheitsmatrix als Startwert für Sensitivitäts-DGL
                    S_b = zeros(n_y);
                    while (i < obj.maxit)
                        % Lösen des Anfangwertproblems
                        sol_y = obj.ode_method(g,tspan,eta,obj.options);
                        y_t_eta = sol_y.y';
                        % Berechnung von F(eta)
                        F = r(eta,y_t_eta(end,:));
                        %
                        Norm_F(i+1) = norm(F);
                        if (norm(F) < obj.StopTol)
                            break; % Abbruchbedingung
                        end
                        
                        % Berechnung der Jacobimatrix 
                        S_t =@(t,S) reshape(g_y(t,deval(sol_y,t))*reshape(S,n_y,n_y),[],1); % Marix-Matrix-Produnkt und zurück in Spaltenvektor 
                        sol_S = obj.ode_method(S_t,tspan,reshape(S_a,[],1),obj.options);
                        S_b(:,:) = reshape(sol_S.y(:,end)',n_y,n_y); % Spaltenvektor zu Matrix
                        
%                         % Berechnung der Jacobimatrix: SCHNELLER, GIBT aber eine Differenz im Bereich 1e-4
%                         for j=1:n_y
%                             odefun = @(time,x) g_y(time,deval(sol_y,time))*x;
%                             [~,S] = obj.ode_method(odefun,tspan,S_a(:,j),obj.options); 
%                             S_b(:,j)=S(end,:)';
%                         end
                        
                        F_jac = r_y_a(eta,y_t_eta(end,:)) + r_y_b(eta,y_t_eta(end,:))*S_b;
                        % Berechnung der Newton-Richtung d
                        d = - F_jac \ F;
                        % Armijo Schrittweiten steuerung
                        while obj.armijo
                            if norm(d) < obj.StopTolArmijo
                                disp('Warning: Step size too small!');
                                t_out = sol_y.x';
                                y_out = y_t_eta;
                                return
                            end
                            [~,y_step] = obj.ode_method(g,tspan,eta+d,obj.options);
                            F_step = r(eta+d,y_step(end,:));
                            if norm(F_step) <= norm(F) + obj.sigma*norm(F_jac*d)
                                break
                            end
                            d = obj.beta*d;
                        end
                        % Neue Werte eta und i
                        eta = eta + d;
                        i = i + 1;
                    end
                case 'FinitDiff'
                    %% Benötigten Ableitungen mit Sensitivitätsdifferentialgleichung lösen
                    I = eye(n_y); % Einheitsmatrix für Einheitsvektoren
                    while (i < obj.maxit)
                        % Lösen des Anfangwertproblems
                        sol_y = obj.ode_method(g,tspan,eta,obj.options);
                        y_t_eta = sol_y.y';  
                        % Berechnung von F(eta)
                        F = r(eta,y_t_eta(end,:));
                        %
                        Norm_F(i+1) = norm(F);
                        if (norm(F) < obj.StopTol)
                            break; % Abbruchbedingung
                        end
                        % Berechnung der Jacobimatrix
                        F_jac = zeros(n_r,n_y);
                        for j=1:n_y
                            y_a_h = eta + obj.h_min * I(:,j); % Erneutes lösen des AWP für ein neues y_a = (eta+h*e_j)
                            sol_y_h = obj.ode_method(g,tspan,y_a_h,obj.options);
                            F_h = r(y_a_h,sol_y_h.y(:,end));
                            F_jac(:,j) = (F_h-F)/obj.h_min; % bilden der Jac-Matrix über finite Differenzen
                        end
                        % Berechnung der Newton-Richtung d 
                        d = - F_jac \ F;
                        % Armijo Schrittweitensteuerung
                        while obj.armijo
                            if norm(d) < obj.StopTolArmijo
                                disp('Warning: Step size too small!');
                                t_out = sol_y.x';
                                y_out = y_t_eta;
                                return
                            end
                            [~,y_step] = obj.ode_method(g,tspan,eta+d,obj.options);
                            F_step = r(eta+d,y_step(end,:));
                            if norm(F_step) <= norm(F) + obj.sigma*norm(F_jac*d)
                                break
                            end
                            d = obj.beta*d;
                        end
                        % Neue Werte eta und i
                        eta = eta + d;
                        i = i + 1;
                    end
                otherwise
                    fprintf(2,'ERROR: Falsches Verfahren für die Berechnung der Jacobimatrix!\n');
                    return;
            end
            if i == obj.maxit
                disp('Warning: Max iterations reached!');
            end
            % Rückgabewerte
            t_out = sol_y.x';
            y_out = y_t_eta;
        end
        
        %% Multiple shooting method
        function [t_out,y_out,eta,i,Norm_F] = Mehrfachschiessverfahren(obj,tspan,eta_0,g,g_y,r,r_y_a,r_y_b)
            % Parameter und Einstellungen
            N = size(eta_0,2);
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
            S_b = zeros(n_y);

            function F_Fjac_func(eta)
                % Zurücksetzen der alten Werte
                t_out = [];
                y_out = [];
                % AWP lösen und aufstellen von F(eta) und der Jacobimatrix
                for j = 1:N
                    % ODE und Sensitivitäts DGL lösen 
                    sol_y = obj.ode_method(g,[t(j),t(j+1)],eta(:,j),obj.options);
                    y_AWP = sol_y.y';
                    t_AWP = sol_y.x';

                    S_t =@(t,S) reshape(g_y(t,deval(sol_y,t))*reshape(S,n_y,n_y),[],1); % Marix-Matrix-Produnkt und zurück in Spaltenvektor 
                    sol_S = obj.ode_method(S_t,[t(j),t(j+1)],reshape(S_a,[],1),obj.options);
                    S_b_AWP(:,:) = reshape(sol_S.y(:,end)',n_y,n_y); % Spaltenvektor zu Matrix

%                     % Berechnung der Jacobimatrix: SCHNELLER, GIBT aber eine Differenz im Bereich 1e-4
%                     for k=1:n_y
%                         odefun = @(time,x) g_y(time,deval(sol_y,time))*x;
%                         [~,S] = obj.ode_method(odefun,[t(j),t(j+1)],S_a(:,k),obj.options); 
%                         S_b_AWP(:,k)=S(end,:)';
%                     end
                    
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
            while (i < obj.maxit)
                % AWP lösen und aufstellen von F(eta) und der Jacobimatrix
                F_Fjac_func(eta);
                % Abbruchbedingung
                Norm_F(i+1) = norm(F);
                if (norm(F) <= obj.StopTol)
                    break; 
                end
                % Berechnung der Newton-Richtung d 
                d = - F_jac \ F;
                d_r = reshape(d,n_y,N);
                % Armijo-Regel zur Schrittweitensteuerung
                t_armijo = 1; % initial t
                F_Fjac_func(eta + t_armijo * d_r);
                while (F > (F + obj.sigma * t_armijo * F_jac * d))
                    t_armijo = obj.beta*t_armijo;
                    F_Fjac_func(eta + t_armijo * d_r);
                end
                % Neue Werte eta und i
                eta = eta + t_armijo*d_r;
                i = i + 1;
            end
        end
    end
end
